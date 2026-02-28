# Feedback Driven Wind Simulations - AREPO
AREPO version that was developed for the purpose of building stellar feedback driven wind simulations. The code was developed for a Master Thesis at the Heidelberg Institute of Theoretical Astrophysics under the supervision of Dr. Dylan Nelson. The original code was developed by Volker Springel et al. and can be found [here](https://gitlab.mpcdf.mpg.de/vrs/arepo). This modified version is distributed under the GNU GPL license.

Click [here](https://github.com/nelson-group/raph_thesis_scripts/blob/main/msc_thesis.pdf) to view a copy of the Master Thesis. 

The code contains several new configuration options to simulate the wind physics. The core wind physics are configured using `INJECT_WITHIN_RADIUS`, which models the dynamics of a starburst-driven wind based on parameterizations of the analytic wind model of Chevalier and Clegg 1985, or [CC85](https://www.nature.com/articles/317044a0). At each time-step in the simulation, mass and energy is deposited in a spherical starburst region with $R_{\rm inject}$. The distribution of mass and energy can be using either:
1) (Default) An even distribution of mass and energy across each cell contained inside a radius of $r \leq R_{\rm inject}$.
2) (`VOLUME_BASED_INJECTION`), where cells are weighted by their volume and normalized by the total volume of all cells where $r \leq R_{\rm inject}$.
The injection process can be controlled by several key parameters, such as injection radius $R_{\rm inject}$, mass load $\beta$, energy load $\alpha$, and burst duration $t_{\rm burst}$. 

Gravitational effects of the disk on the wind flow is represented in the form of a [Miyamoto-Nagai potential](https://ui.adsabs.harvard.edu/abs/1975PASJ...27..533M/abstract), which can be controlled using by configuring stellar disk mass `M_stars`, stellar scale radius `R_stars`, and stellar scale height `z_stars`

Additionally, the code also implements a simple metallic radiative cooling scheme that reads off of publically available lookup tables of [Wiersma et al. 2009](https://local.strw.leidenuniv.nl/WSS08/). Two options for mettalic cooling are included: (i) `CIE_PIE_COOLING`, which is dependent on both density and temperature, and (ii) `CIE_COOLING`, which is solely temperature dependent. 

Avenues for future development include implementing more complex stellar feedback routines (For example, a cluster feedback model similar to the one implemented by [Schneider et al. 2018](https://arxiv.org/abs/1803.01005)), as a well a more extensive radiative cooling model. 