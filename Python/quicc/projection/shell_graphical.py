# -*- coding: utf-8 -*-
""" Implementation of shell_graphical.
author: Nicolò Lardelli
data: 07.02.18
"""
import h5py
from quicc.projection import spherical, shell
from quicc.geometry.spherical  import shell_radius as geo
import numpy as np
from numpy.polynomial import chebyshev as cheb
from numpy.polynomial import legendre as leg
from matplotlib import pyplot as pp
from matplotlib import rc, rcParams
from matplotlib.colors import LogNorm, Normalize

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
pp.rc('text.latex', preamble=r'\usepackage{bm}')

def rank_1_matrix(a, b):
    # Input:
    # a: column vector
    # b: row vector

    a = np.reshape(a, (-1, 1))
    b = np.reshape(b, (-1, 1))

    return np.kron(a, b.T)


class ShellPlotter:

    def __init__(self, filename):

        # open file and object
        self.fopen = h5py.File(filename, 'r')

        # determine type
        self.file_type = self.fopen.attrs['type']

        # retrieve parameters for visualization
        self.ro = self.fopen['/physical/ro'].value
        self.rratio = self.fopen['/physical/rratio'].value
        self.E = self.fopen['/physical/ekman'].value

        # retrieve the runtime and the frequency
        self.time = self.fopen['/run/time'].value
        try:
            self.fN = self.fopen['/physical/omega'].value
        except:
            # assume there are no time depenencies
            self.fN = 0

        # define the phi_0 variable (meriodional point used to plot meridional sections)
        self.phi_0 = self.fN * self.time

        # retrieve the spectral resolution
        self.nN = self.fopen['truncation/spectral/dim1D'].value + 1
        self.nL = self.fopen['truncation/spectral/dim2D'].value + 1
        self.nM = self.fopen['truncation/spectral/dim3D'].value + 1

        # produce the mapping
        self.a, self.b = geo.linear_r2x(self.ro, self.rratio)
        self.ri = self.ro * self.rratio

        # define title dictionary
        self.title_dict = {'simple': r'$\bm{u}_', 'curl': r'$\left(\bm{\nabla}\times\bm{u}\right)_'}


    def plot(self, section_type, m = 'all', type = 'simple'):

        # prepare the modes
        if m == 'all':
            modes = np.arange(self.nM)
        else:
            modes = np.array(m)

        # store the field type
        self.vector_field_type = type

        if section_type == 'meridional':
            self.plot_meridional(modes = modes)
        elif section_type == "equatorial":
            self.plot_equatorial(modes = modes)
        elif section_type == "boundaries":
            self.plot_boundaries(modes = modes)
        else:
            raise RuntimeError('Unkown section type '+ section_type)

    def loop_over(self, *args, **kwargs ):

        # loop_over takes care of the looping of the various modes and distinguishes between SLFm and SLFl

        if self.file_type == b'SLFl':

            idx = 0
            for l in range(self.nL):

                for m in range(min(l + 1, self.nM)):

                    if m in kwargs['modes']:
                        self.evaluate_mode(l, m, idx, *args, **kwargs)

                    idx += 1

        elif self.file_type == b'SLFm':

            idx = 0
            for m in range(self.nM):

                for l in range(m, self.nL):

                    if m in kwargs['modes']:
                        self.evaluate_mode(l, m, idx, *args, **kwargs)

                    idx += 1
        else:
            raise RuntimeError('Unknown  file type ' + self.file_type)

    def plot(self,m='all', **kwargs):


        # prepare the modes
        if m == 'all':
            modes = np.arange(self.nM)
        else:
            modes = np.array(m)

        kwargs['modes']=modes
        # produce the grid for plotting
        xx_r, ww = cheb.chebgauss(2 * self.nN)
        rr = self.a * xx_r + self.b

        # produce grid for bulk of the flow only
        delta = self.E ** .5

        # first if decision block over the radial grid
        if kwargs['mode']=='boundaries':


            rr = rr[(rr < self.ri + 10 * delta) | (rr > self.ro - 10 * delta)]
            xx_bulk = (rr - self.b) / self.a

        elif  kwargs['mode']=='line':

            eta = self.rratio
            ymin = eta/(1-eta)
            xmin = 0.
            theta_crit = np.arccos(np.abs(self.fN)/2)
            h = ((1+eta)/(1-eta))**.5
            xmax = xmin + np.sin(theta_crit)*h
            ymax = ymin + np.cos(theta_crit)*h

            xx = np.linspace(xmin, xmax, 4*self.nN )
            yy = np.linspace(ymin, ymax, 4*self.nN )

            rr = (xx**2+yy**2)**.5
            ttheta = np.arctan2(xx,yy)

            # select the interior of the flow
            idx = (rr >= self.ri + 10 * delta) & (rr <= self.ro - 10 * delta)
            rr = rr[idx]

            xx_bulk = (rr - self.b) / self.a
            ttheta = ttheta[idx]
        else:

            #delta = 0.
            rr = rr[(rr >= self.ri + 10 * delta) & (rr <= self.ro - 10 * delta)]
            xx_bulk = (rr - self.b) / self.a


        # second if decision block over the secundary grid
        kwargs['pphi']=None
        if kwargs['mode']=='equatorial':

            # generate the azimuthal grid
            pphi = np.linspace(0, 2 * np.pi, 2 * self.nL + 1) + self.phi_0
            kwargs['pphi'] = pphi
        elif kwargs['mode']=='line':

            self.xx_th = np.cos(ttheta)
            pass
        else:

            # generate the meridional grid
            xx_th, ww = leg.leggauss(2*self.nL)
            ttheta = np.arccos(xx_th)
            ttheta = np.concatenate(([np.pi], ttheta, [0]))
            self.xx_th = np.cos(ttheta)


        # third decision block, for the x/y grid generation
        if kwargs['mode']=='equatorial':

            RR, PP = np.meshgrid(rr, pphi)
            XX = np.cos(PP) * RR
            YY = np.sin(PP) * RR

        elif kwargs['mode']=='line':
            XX=rr
            pass
        elif kwargs['mode'] =='boundaries':

            RR, TT = np.meshgrid(rr, ttheta)
            XX=TT
            YY=RR
        else:

            RR, TT = np.meshgrid(rr, ttheta)
            YY = np.cos(TT) * RR
            XX = np.sin(TT) * RR

        # produce the mapping for chebyshev polynomials
        self.Tn_eval = shell.proj_radial(self.nN, self.a, self.b, xx_bulk)  # evaluate chebyshev simple
        self.dTndr_eval = shell.proj_dradial_dr(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r d/dr(r Tn)
        self.Tn_r_eval = shell.proj_radial_r(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r Tn

        # produce the mapping for the tri-curl part
        self.Tn_r2_eval = shell.proj_radial_r2(self.nN, self.a, self.b, xx_bulk) # evaluate 1/r**2 Tn
        self.d2Tndr2_eval = shell.proj_lapl(self.nN, self.a, self.b, xx_bulk) # evaluate 1/r**2dr r**2 dr

        if kwargs['mode']=='boundaries':

            # prepare the subplot
            fig, ax = pp.subplots(2, 3, sharey=False, sharex=True, figsize=(10, 2))
            ax1 = ax[:,0]
            ax2 = ax[:,1]
            ax3 = ax[:,2]
        else:
            # prepare the subplot
            fig, (ax1, ax2, ax3) = pp.subplots(1, 3, sharey=True, sharex=True, figsize=(10, 4))

        if kwargs['mode'] == 'line':
            x_arg = rr
        else:
            x_arg = (XX, YY)

        if kwargs['type']=='simple' or kwargs['type']=='curl':

            prefix = self.title_dict[kwargs['type']]

            # initialize the fields
            U_phi = np.zeros_like(XX, dtype=complex)
            U_r = np.zeros_like(XX, dtype=complex)
            U_theta = np.zeros_like(XX, dtype=complex)

            # compute (either velocity or vorticity
            self.loop_over(U_r, U_theta, U_phi, **kwargs)


            # plot r field component
            self.plot_field(fig, ax1, x_arg, U_r, colormap='dual', title=prefix+'{r}$', **kwargs)
            if kwargs['mode']!='boundaries' and kwargs['mode']!='line':
                ax1.set_aspect('equal')

            # plot theta field component
            self.plot_field(fig, ax2, x_arg, U_theta, colormap='dual', title=prefix+'{\Theta}$', **kwargs)

            # plot phi field component
            self.plot_field(fig, ax3, x_arg, U_phi, colormap='dual', title=prefix+'{\phi}$', **kwargs)

        else:

            # initialize the fields
            U_phi = np.zeros_like(XX, dtype=complex)
            U_r = np.zeros_like(XX, dtype=complex)
            U_theta = np.zeros_like(XX, dtype=complex)

            # compute velocity
            kwargs['type']='simple'
            self.loop_over(U_r, U_theta, U_phi,  **kwargs)

            # initialize the fields
            Omega_phi = np.zeros_like(XX, dtype=complex)
            Omega_r = np.zeros_like(XX, dtype=complex)
            Omega_theta = np.zeros_like(XX, dtype=complex)

            # compute vorticity
            kwargs['type']='curl'
            self.loop_over(Omega_r, Omega_theta, Omega_phi,  **kwargs)

            # take only real parts
            U_r = np.real(U_r)
            U_theta = np.real(U_theta)
            U_phi = np.real(U_phi)
            Omega_r = np.real(Omega_r)
            Omega_theta = np.real(Omega_theta)
            Omega_phi = np.real(Omega_phi)

            Energy = U_r**2 + U_theta**2 + U_phi**2
            Energy *= 0.5

            Enstrophy = Omega_r**2 +  Omega_theta**2 + Omega_phi**2
            Enstrophy *= 0.5

            Helicity = U_r * Omega_r + U_theta * Omega_theta + U_phi * Omega_phi

            # plot Energy
            self.plot_field(fig, ax1, x_arg, Energy, colormap = 'log', title=r'$\frac{1}{2}\left|\bm{u}\right|^2$', **kwargs)
            if kwargs['mode']!='boundaries' and kwargs['mode']!='line':
                ax1.set_aspect('equal')
            # plot Enstrophy
            self.plot_field(fig, ax2, x_arg, Enstrophy, colormap='log', title=r'$\frac{1}{2}\left|\bm{\nabla}\times\bm{u}\right|^2$', **kwargs)

            # plot Helicity
            self.plot_field(fig, ax3, x_arg, Helicity, colormap='dual', title=r'$\bm{u}\cdot\bm{\nabla}\times\bm{u}$', **kwargs)

    """
    def plot_equatorial(self, **kwargs):

        # produce the grid for plotting
        xx_r, ww = cheb.chebgauss(4 * self.nN)
        rr = self.a * xx_r + self.b

        # produce grid for bulk of the flow only
        # delta = self.E ** .5
        delta = 0.
        rr = rr[(rr >= self.ri + 10 * delta) & (rr <= self.ro - 10 * delta)]
        xx_bulk = (rr - self.b) / self.a

        # generate the grid for the
        pphi = np.linspace(0, 2*np.pi, 2 * self.nL+1) + self.phi_0
        RR, PP = np.meshgrid(rr, pphi)
        XX = np.cos(PP) * RR
        YY = np.sin(PP) * RR

        # produce the mapping for chebyshev polynomials
        self.Tn_eval = shell.proj_radial(self.nN, self.a, self.b, xx_bulk)  # evaluate chebyshev simple
        self.dTndr_eval = shell.proj_dradial_dr(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r d/dr(r Tn)
        self.Tn_r_eval = shell.proj_radial_r(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r Tn

        # produce the mapping for the tri-curl part
        self.Tn_r2_eval = shell.proj_radial_r2(self.nN, self.a, self.b, xx_bulk) # evaluate 1/r**2 Tn
        self.d2Tndr2_eval = shell.proj_lapl(self.nN, self.a, self.b, xx_bulk) # evaluate 1/r**2dr r**2 dr

        # initialize the fields
        U_phi = np.zeros_like(XX, dtype=complex)
        U_r = np.zeros_like(XX, dtype=complex)
        U_theta = np.zeros_like(XX, dtype=complex)

        # loop over
        self.loop_over(U_r, U_theta, U_phi, pphi = pphi, **kwargs)

        if kwargs['mode']=='boundaries':

            # prepare the subplot
            fig, ax = pp.subplots(2, 3, sharey=False, sharex=True, figsize=(10, 2))
            ax1 = ax[:,0]
            ax2 = ax[:,1]
            ax3 = ax[:,2]
        else:
            # prepare the subplot
            fig, (ax1, ax2, ax3) = pp.subplots(1, 3, sharey=True, sharex=True, figsize=(10, 4))

        prefix = self.title_dict[self.vector_field_type]

        # plot Ur
        im1 = ax1.contourf(XX, YY, np.real(U_r), 100, cmap = pp.get_cmap('coolwarm'), vmin=np.nanmin(np.real(U_r)),
                           vmax=np.nanmax(np.real(U_r)))
        ax1.set_aspect('equal')
        ax1.set_title(prefix+'{r}$')
        pp.colorbar(im1, ax = ax1, orientation = 'horizontal', ticks=[np.nanmin(np.real(U_r)), np.nanmax(np.real(U_r))]
                    , shrink=0.8)

        # plot Utheta
        im2 = ax2.contourf(XX, YY, np.real(U_theta), 100, cmap = pp.get_cmap('coolwarm'), vmin=np.nanmin(np.real(U_theta)),
                           vmax=np.nanmax(np.real(U_theta)))
        ax2.set_title(prefix+'{\theta}$')
        pp.colorbar(im2, ax = ax2, orientation = 'horizontal',
                    ticks=[np.nanmin(np.real(U_theta)), np.nanmax(np.real(U_theta))], shrink=0.8)

        # plot Uphi
        im3 = ax3.contourf(XX, YY, np.real(U_phi), 100, cmap = pp.get_cmap('coolwarm'), vmin=np.nanmin(np.real(U_phi)),
                           vmax=np.nanmax(np.real(U_phi)))
        ax3.set_title(prefix+'{\phi}$')
        pp.colorbar(im3, ax = ax3, orientation = 'horizontal',
                    ticks=[np.nanmin(np.real(U_phi)), np.nanmax(np.real(U_phi))], shrink=0.8)
        print(np.nanmin(np.real(U_phi)), np.nanmax(np.real(U_phi)))
    """
    """

    def plot_boundaries(self, **kwargs):

        # produce the grid for plotting
        xx_r, ww = cheb.chebgauss(2 * self.nN)
        rr = self.a * xx_r + self.b

        xx_th, ww = leg.leggauss(2 * self.nL)
        ttheta = np.arccos(xx_th)

        # produce grid for the boundary
        delta = self.E ** .5
        rr = rr[(rr < self.ri + 10 * delta) | (rr > self.ro - 10 * delta)]
        xx_bulk = (rr - self.b) / self.a
        ttheta = np.concatenate(([np.pi], ttheta, [0]))
        self.xx_th = np.cos(ttheta)

        RR, TT = np.meshgrid(rr, ttheta)

        # produce the mapping for chebyshev polynomials
        self.Tn_eval = shell.proj_radial(self.nN, self.a, self.b, xx_bulk)  # evaluate chebyshev simple
        self.dTndr_eval = shell.proj_dradial_dr(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r d/dr(r Tn)
        self.Tn_r_eval = shell.proj_radial_r(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r Tn

        # produce the mapping for the tri-curl part
        self.Tn_r2_eval = shell.proj_radial_r2(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r**2 Tn
        self.d2Tndr2_eval = shell.proj_lapl(self.nN, self.a, self.b, xx_bulk)  # evaluate 1/r**2dr r**2 dr

        # initialize the fields
        U_phi = np.zeros_like(RR, dtype=complex)
        U_r = np.zeros_like(RR, dtype=complex)
        U_theta = np.zeros_like(RR, dtype=complex)

        # compute U_r, U_theta, U_phi
        self.loop_over(U_r, U_theta, U_phi, **kwargs)

        # prepare the subplot
        fig, ax = pp.subplots(2,3, sharey=False, sharex=True, figsize = (10,2))
        # find the separation between inner and outer boundary
        idx_boundary = int(TT.shape[1]/2)

        fig.subplots_adjust(hspace=0)

        # define dictionary to plot
        prefix = self.title_dict[self.vector_field_type]
        # plot Ur

        im00 = ax[0, 0].contourf(TT[:, 0:idx_boundary], RR[:, 0:idx_boundary], np.real(U_r[:, 0:idx_boundary]), 100,
                                 cmap = pp.get_cmap('coolwarm'), vmin=min, vmax=max)
        ax[0, 0].set_aspect(15.0)
        im10 = ax[1, 0].contourf(TT[:, idx_boundary + 1:], RR[:, idx_boundary + 1:], np.real(U_r[:, idx_boundary + 1:]),
                                 100, cmap=pp.get_cmap('coolwarm'), vmin=min, vmax=max)
        ax[1, 0].set_aspect(15.0)
        ax[0, 0].set_title(prefix+r'{r}$')
        cb1 = fig.colorbar(im10, orientation='horizontal', ax=[ax[0, 0], ax[1, 0]], ticks=[min, max], shrink=0.8)

        # plot Utheta
        min = np.nanmin(np.real(U_theta))
        max = np.nanmax(np.real(U_theta))
        im01 = ax[0,1].contourf(TT[:, 0:idx_boundary], RR[:, 0:idx_boundary], np.real(U_theta[:, 0:idx_boundary]), 100,
                                cmap=pp.get_cmap('coolwarm'), vmin=min, vmax=max)
        ax[0, 1].set_aspect(15.0)
        im11 = ax[1,1].contourf(TT[:, idx_boundary + 1:], RR[:, idx_boundary + 1:],
                                np.real(U_theta[:, idx_boundary + 1:]), 100, cmap=pp.get_cmap('coolwarm'),
                                vmin=min, vmax=max)
        ax[1, 1].set_aspect(15.0)
        ax[0,1].set_title(prefix+r'{\theta}$')
        #print(np.real(U_theta).nanmin(), np.real(U_theta).nanmax())
        cb2 = fig.colorbar(im11, orientation='horizontal', ax = [ax[0, 1], ax[1, 1]], ticks = [min, max], shrink=0.8)

        # plot Uphi
        min = np.nanmin(np.real(U_phi))
        max = np.nanmax(np.real(U_phi))
        im02 = ax[0,2].contourf(TT[:, 0:idx_boundary], RR[:, 0:idx_boundary], np.real(U_phi[:, 0:idx_boundary]), 100,
                               cmap=pp.get_cmap('coolwarm'), vmin=min, vmax=max)
        ax[0, 2].set_aspect(15.0)
        im12 = ax[1,2].contourf(TT[:, idx_boundary + 1:], RR[:, idx_boundary + 1:], np.real(U_phi[:, idx_boundary + 1:]),
                               100, cmap=pp.get_cmap('coolwarm'), vmin=min, vmax=max)
        ax[1, 2].set_aspect(15.0)
        ax[0, 2].set_title(prefix+r'{\phi}$')
        cb3 = fig.colorbar(im12, orientation='horizontal', ax = [ax[0, 2], ax[1, 2]], ticks =[min, max], shrink=0.8)

    """
    def plot_field(self, fig, ax, *args, **kwargs):
        if kwargs['mode']=='line':
            rr = args[0]
            ff = args[1]

            ff = np.real(ff)
            min = np.nanmin(ff)
            max = np.nanmax(ff)
        else:
            XX = args[0]
            YY = args[1]
            ZZ = args[2]
            ZZ = np.real(ZZ)
            min = np.nanmin(ZZ)
            max = np.nanmax(ZZ)

        if kwargs.get('colormap','dual') == 'log':
            map = pp.get_cmap('hot')
            norm = LogNorm(min, max)

        else:
            map = pp.get_cmap('coolwarm')
            norm = Normalize(min, max, clip=True)

        if kwargs['mode']=='boundaries':
            idx_boundary = idx_boundary = int(XX.shape[1]/2)

            ax[0].contourf(XX[:, 0:idx_boundary], YY[:, 0:idx_boundary], ZZ[:, 0:idx_boundary], 100, cmap=map, norm = norm)
            ax[0].set_aspect(15.0)
            im = ax[1].contourf(XX[:, idx_boundary + 1:], YY[:, idx_boundary + 1:], ZZ[:, idx_boundary + 1:], 100, cmap=map, norm = norm)
            ax[1].set_aspect(15.0)
            ax[0].set_title(kwargs['title'])
            #cb1 = fig.colorbar(im10, orientation='horizontal', ax=[ax[0], ax[1]], ticks=[min, max], shrink=0.8)
        elif kwargs['mode']=='line':
            ax.plot(rr, ff )
        else:
            im = ax.contourf(XX, YY, ZZ, 100, cmap=map, norm = norm)
            ax.set_title(kwargs['title'])
            if kwargs.get('colormap', 'dual') == 'log':
                fig.colorbar(im, orientation='horizontal', ax=ax, shrink=1.0)
            else:
                fig.colorbar(im, orientation='horizontal', ax=ax, shrink=1.0)



    def evaluate_mode(self, l, m, idx, *args, **kwargs):

        Field_r = args[0]
        Field_theta = args[1]
        Field_phi = args[2]

        # retrieve the scalar fields
        modeP = self.fopen['velocity/velocity_pol'].value[idx, :, 0] + self.fopen['velocity/velocity_pol'].value[idx, :, 1] * 1j
        modeT = self.fopen['velocity/velocity_tor'].value[idx, :, 0] + self.fopen['velocity/velocity_tor'].value[idx, :, 1] * 1j

        # initialize radial parts
        # the radial parts can be interchanged if we evaluate the vorticity field or the velocity field
        # we need to invert toroidal with poloidal however
        if kwargs['type'] == 'simple':
            # procedure for the velocity field
            rad_part_ur = l * (l + 1) * np.matmul(self.Tn_r_eval, modeP)
            rad_part_pol = np.matmul(self.dTndr_eval, modeP)
            rad_part_tor = np.matmul(self.Tn_eval, modeT)

        elif kwargs['type'] == 'curl':
            # procedure for the vorticity field
            rad_part_ur = l * (l + 1) * np.matmul(self.Tn_r_eval, modeT)
            rad_part_pol = np.matmul(self.dTndr_eval, modeT)
            rad_part_tor = -( np.matmul(self.d2Tndr2_eval, modeP)+ l * (l+1) *np.matmul(self.Tn_r2_eval, modeP) )

        else:
            raise RuntimeError('Unknown vector field type '+self.vector_field_type)
            pass

        factor = 2.0

        if kwargs['mode']=='meridional':

            # prepare arrays
            eimp = np.exp(1j * m * self.phi_0)*factor

            # initialize the theta tranforms
            Ylm = spherical.lplm(self.nL, l, m, self.xx_th)
            dYlm = spherical.dplm(self.nL, l, m, self.xx_th)
            Ylm_sin = spherical.lplm_sin(self.nL, l, m, self.xx_th) * m * 1j

            # update the fields poloidal parts
            Field_r += rank_1_matrix(Ylm, rad_part_ur) * eimp
            Field_theta += rank_1_matrix(dYlm, rad_part_pol) * eimp
            Field_phi += rank_1_matrix(Ylm_sin, rad_part_pol) * eimp

            # update the fields toroidal parts
            Field_theta += rank_1_matrix(Ylm_sin, rad_part_tor) * eimp
            Field_phi -= rank_1_matrix(dYlm, rad_part_tor) * eimp

        elif kwargs['mode']=='line':

            # prepare arrays
            eimp = np.exp(1j * m * self.phi_0)*factor

            # initialize the theta tranforms
            Ylm = spherical.lplm(self.nL, l, m, self.xx_th)
            dYlm = spherical.dplm(self.nL, l, m, self.xx_th)
            Ylm_sin = spherical.lplm_sin(self.nL, l, m, self.xx_th) * m * 1j

            # update the fields poloidal parts
            Field_r += Ylm* rad_part_ur * eimp
            Field_theta += dYlm* rad_part_pol * eimp
            Field_phi += Ylm_sin* rad_part_pol * eimp

            # update the fields toroidal parts
            Field_theta += Ylm_sin* rad_part_tor * eimp
            Field_phi -= dYlm* rad_part_tor * eimp
        else:

            # prepare arrays
            eimp = np.exp(1j * m * kwargs['pphi'])*factor

            # initialize the theta tranforms
            Ylm = spherical.lplm(self.nL, l, m, np.array([0.]))
            dYlm = spherical.dplm(self.nL, l, m,  np.array([0.]))
            Ylm_sin = spherical.lplm_sin(self.nL, l, m,  np.array([0.])) * m * 1j

            # update the fields poloidal parts
            Field_r += rank_1_matrix(eimp, rad_part_ur) * Ylm
            Field_theta += rank_1_matrix(eimp, rad_part_pol) * dYlm
            Field_phi += rank_1_matrix(eimp, rad_part_pol) * Ylm_sin

            # update the fields toroidal parts
            Field_theta += rank_1_matrix(eimp, rad_part_tor) * Ylm_sin
            Field_phi -= rank_1_matrix(eimp, rad_part_tor) * dYlm