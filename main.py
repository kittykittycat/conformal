import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, ifft2, fftfreq

from tkinter import *
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

res = 48
zoom = 1.1
num_zeros = 4
num_poles = 4
r = 4

gridx, gridy = np.mgrid[-1:1:res*1j, -1:1:res*1j]
grid = gridx + 1j*gridy
#grid = np.exp(grid)

rng = np.random.default_rng()
zeros = np.linspace(0, 1, num=num_zeros, endpoint=False) + rng.random()
poles = np.linspace(0, 1, num=num_poles, endpoint=False) + rng.random()
zeros = r*np.exp(2j*np.pi*zeros)
poles = r*np.exp(2j*np.pi*poles)

def trafo(x):
    result = 1
    for z in zeros:
        result *= x-z
    for p in poles:
        result /= x-p
    m = np.abs(result).max()
    return result/m

image = trafo(grid)
h_x = image.real
h_y = image.imag
v_x = image.T.real
v_y = image.T.imag




root = Tk()
root.title("conformal function plotter")

frm_toplevel = ttk.Frame(root, padding=10)
frm_toplevel.grid()

fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
pl1 = ax.plot(h_x, h_y, "b", dashes=[6, 2], linewidth=0.6)
pl2 = ax.plot(v_x, v_y, "r", dashes=[2, 1], linewidth=0.6)
pl3, = ax.plot(zeros.real, zeros.imag, "k.")
pl4, = ax.plot(poles.real, poles.imag, "g.")
ax.set_xlim([-zoom,zoom])
ax.set_ylim([-zoom,zoom])
canvas = FigureCanvasTkAgg(fig, master=frm_toplevel)
canvas.draw()


def reroll():
    global zeros
    global poles
    global image
    zeros = np.linspace(0, 1, num=num_zeros, endpoint=False) + rng.random()
    zeros = r*np.exp(2j*np.pi*zeros)
    poles = np.linspace(0, 1, num=num_poles, endpoint=False) + rng.random()
    poles = r*np.exp(2j*np.pi*poles)
    image = trafo(grid)
    pl3.set_data(zeros.real, zeros.imag)
    pl4.set_data(poles.real, poles.imag)

def update(val):
    val = float(val)
    val = val**(8/(num_zeros**2+num_poles**2+r**2))

    # convex combo of preimage and image
    img = (1-val)*grid + val*image

    # scale such that we use the [-1,1]x[-i,i] square
    img -= img.real.min() + 1j*img.imag.min()
    img *= 2/max( img.real.max(), img.imag.max() )
    img -= 1+1j

    # the grid lines
    hxnew = img.real
    hynew = img.imag
    vxnew = img.T.real
    vynew = img.T.imag

    # update plot data
    for k, line in enumerate(pl1):
        line.set_data(hxnew[k],
                      hynew[k])

    for k, line in enumerate(pl2):
        line.set_data(vxnew[k],
                      vynew[k])
    canvas.draw()

# pack_toolbar=False will make it easier to use a layout manager later on.
toolbar = NavigationToolbar2Tk(canvas, frm_toplevel, pack_toolbar=False)
toolbar.update()
toolbar.pack(side=BOTTOM, fill=X)
canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)

slider = Scale(frm_toplevel,
               from_=0.0,
               to=1.0,
               orient=HORIZONTAL,
               length=500,
               command=update,
               showvalue=False,
               #digits=2,
               resolution=1/500)
slider.set(0)
slider.pack(side=BOTTOM)
update(0)

btn_reroll = ttk.Button(frm_toplevel, text="reroll", command=reroll)
btn_reroll.pack(side=BOTTOM)

root.mainloop()