


# import matplotlib.pyplot as plt
# import numpy as np

# w = 3
# Y, X = np.mgrid[-w:w:100j, -w:w:100j]
# U = -1 - X**2 + Y
# V = 1 + X - Y**2
# speed = np.sqrt(U**2 + V**2)

# fig, axs = plt.subplots(3, 2, figsize=(7, 9), height_ratios=[1, 1, 2])
# axs = axs.flat

# #  Varying density along a streamline
# axs[0].streamplot(X, Y, U, V, density=[0.5, 1])
# axs[0].set_title('Varying Density')

# # Varying color along a streamline
# strm = axs[1].streamplot(X, Y, U, V, color=U, linewidth=2, cmap='autumn')
# fig.colorbar(strm.lines)
# axs[1].set_title('Varying Color')

# #  Varying line width along a streamline
# lw = 5*speed / speed.max()
# axs[2].streamplot(X, Y, U, V, density=0.6, color='k', linewidth=lw)
# axs[2].set_title('Varying Line Width')

# # Controlling the starting points of the streamlines
# seed_points = np.array([[-2, -1, 0, 1, 2, -1], [-2, -1,  0, 1, 2, 2]])

# strm = axs[3].streamplot(X, Y, U, V, color=U, linewidth=2,
#                          cmap='autumn', start_points=seed_points.T)
# fig.colorbar(strm.lines)
# axs[3].set_title('Controlling Starting Points')

# # Displaying the starting points with blue symbols.
# axs[3].plot(seed_points[0], seed_points[1], 'bo')
# axs[3].set(xlim=(-w, w), ylim=(-w, w))

# # Create a mask
# mask = np.zeros(U.shape, dtype=bool)
# mask[40:60, 40:60] = True
# U[:20, :20] = np.nan
# U = np.ma.array(U, mask=mask)

# axs[4].streamplot(X, Y, U, V, color='r')
# axs[4].set_title('Streamplot with Masking')

# axs[4].imshow(~mask, extent=(-w, w, -w, w), alpha=0.5, cmap='gray',
#               aspect='auto')
# axs[4].set_aspect('equal')

# axs[5].streamplot(X, Y, U, V, broken_streamlines=False)
# axs[5].set_title('Streamplot with unbroken streamlines')

# plt.tight_layout()
# plt.show()


import numpy as np
a = np.full((4,4), 2)
b = np.full((4,4), 3)
b[2,3] = 6
# a[:-1,-1] = 100
c = a[1:-1,:]*b[1:-1,:]
print(c)