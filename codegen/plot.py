from matplotlib import pyplot as plt
import numpy as np
import sys
import os

eps = 2.2204e-16

def main():
    if (len(sys.argv) < 5):
        print('usage: {} nt Nx Ny Nz'.format(sys.argv[0]))
        sys.exit(1)
    nt = int(sys.argv[1])
    Nx = int(sys.argv[2])
    Ny = int(sys.argv[3])
    Nz = int(sys.argv[4])

    cx = int(Nx/2)
    cy = int(Ny/2)
    cz = int(Nz/2)

    filename = 'output.csv'

    fig, axs = plt.subplots(1, 3)
    
    chunk_count = Nx * Ny * Nz
    chunk_size = 8 * chunk_count

    os.makedirs('img', exist_ok=True)

    with open(filename, 'rb') as f:
        for t in range(nt):
            chunk = f.read(chunk_size)
            arr = np.copy(np.frombuffer(chunk, dtype=np.float64, count=chunk_count))
            Ez = arr.reshape((Nx, Ny, Nz))

            Ez[Ez == 0] = eps
            E = 10*np.log10(np.abs(Ez))

            E[np.isinf(E)] = 10*np.log10(eps)

            I0 = E[:, :, cz]
            I1 = E[:, cy, :]
            I2 = E[cx, :, :]

            if t == 0:
                im0 = axs[0].imshow(I0, cmap='jet', vmin=10 * np.log10(eps), vmax=0)
                im1 = axs[1].imshow(I1, cmap='jet', vmin=10 * np.log10(eps), vmax=0)
                im2 = axs[2].imshow(I2, cmap='jet', vmin=10 * np.log10(eps), vmax=0)
                plt.show(block=False)

            im0.set_data(I0)
            im1.set_data(I1)
            im2.set_data(I2)
            fig.canvas.draw()
            plt.suptitle('t = {}'.format(t))

            #plt.pause(0.01)
            print('t={}/{}'.format(t, nt-1), end='\r')
            plt.savefig('img/i.{:05}.png'.format(t), dpi=300)
        print()

    # try to create mp4 using ffmpeg if it installed
    try:
        os.system('ffmpeg -y -framerate 10 -pattern_type glob -i "img/*.png" -c:v libx264 -pix_fmt yuv420p sim.mp4')
    except:
        print('error: failed to create mp4 from img/ files, perhaps an issue with ffmpeg.')

if __name__ == '__main__':
    main()
