import numpy as np
import matplotlib.pyplot as plt

from nicegui import events, ui

def show(event: events.ValueChangeEventArguments):
    name = type(event.sender).__name__
    ui.notify(f'{name}: {event.value}')

class renderer:
    def __init__(self):
        self.res = 48
        self.num_zeros = 4
        self.num_poles = 3
        self.radius = 4

        gridx, gridy = np.mgrid[-1:1:self.res*1j, -1:1:self.res*1j]
        self.grid = gridx + 1j*gridy

        self.hx = gridx
        self.hy = gridy
        self.vx = gridx.T
        self.vy = gridy.T

        self.set_trafo()

    def set_trafo(self, num_zeros=None, num_poles=None):
        self.num_zeros = num_zeros if num_zeros is not None else self.num_zeros
        self.num_poles = num_poles if num_poles is not None else self.num_poles
        zeros = np.linspace(0, 1, num=self.num_zeros, endpoint=False)
        poles = np.linspace(0, 1, num=self.num_poles, endpoint=False)
        zeros = self.radius * np.exp(2j * np.pi * zeros)
        poles = self.radius * np.exp(2j * np.pi * poles)

        def trafo(x):
            result = 1
            for z in zeros:
                result *= x-z
            for p in poles:
                result /= x-p
            m = np.abs(result).max()
            return result/m

        self.trafo = trafo

    def update(self, val):
        val = float(val)
        val = val ** (8 / (self.num_zeros ** 2 + self.num_poles ** 2 + self.radius ** 2))

        # convex combo of preimage and image
        self.image = self.trafo(self.grid)
        img = (1 - val) * self.grid + val * self.image

        # scale such that we use the [-1,1]x[-i,i] square
        img -= img.real.min() + 1j * img.imag.min()
        img *= 2 / max(img.real.max(), img.imag.max())
        img -= 1 + 1j

        # the grid lines
        self.hx = img.real
        self.hy = img.imag
        self.vx = img.T.real
        self.vy = img.T.imag

class rendererGUI:
    slider_res = 500

    def __init__(self, renderer):
        self.rend = renderer
        self.alpha = 0
        self.init_plot()
        ui.run()

    def init_plot(self):
        self.main_plot = ui.matplotlib(figsize=(6, 6)).figure
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
        ax = self.main_plot.gca()
        obj = self.rend
        self.pl1 = ax.plot(obj.hx, obj.hy, "b", dashes=[6, 2], linewidth=0.6)
        self.pl2 = ax.plot(obj.vx, obj.vy, "r", dashes=[2, 1], linewidth=0.6)

        main_slider = ui.slider(min=0, max=self.slider_res, on_change=self.update_plot)
        main_slider.bind_value(self, "alpha")
        slice_pos = ui.label()
        slice_pos.bind_text_from(main_slider, 'value')

    def update_plot(self):
        self.rend.update(self.alpha/self.slider_res)
        with self.main_plot:
            for k, line in enumerate(self.pl1):
                line.set_data(self.rend.hx[k],
                              self.rend.hy[k])

            for k, line in enumerate(self.pl2):
                line.set_data(self.rend.vx[k],
                              self.rend.vy[k])


rend = renderer()
rendGUI = rendererGUI(rend)