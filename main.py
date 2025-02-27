import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from copy import copy

from nicegui import events, ui

functions = {"²" : np.square,
             "sin" : np.sin,
             "cos" : np.cos,
             "tan" : np.tan,
             "arcsin" : np.asin,
             "arccos" : np.acos,
             "arctan" : np.atan,
             "sinh" : np.sinh,
             "cosh" : np.cosh,
             "tanh" : np.tanh,
             "arcsinh" : np.asinh,
             "arccosh" : np.acosh,
             "arctanh" : np.atanh,
             "exp" : np.exp,
             "ln" : np.log,
             "sinc" : sp.sinc,
             "1/" : np.reciprocal,
             "erf" : sp.erf,
             "Γ" : sp.gamma,
             "digamma" : sp.psi,
             "W" : sp.lambertw,
             "exp1" : sp.exp1,
             "Wright Omega" : sp.wrightomega,
             # "" : sp.,
             # "" : sp.,
             "√" : np.sqrt}

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
        #val = val ** (8 / (self.num_zeros ** 2 + self.num_poles ** 2 + self.radius ** 2))

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
        self.available_trafos = copy(functions)
        self.available_trafos["rational"] = copy(rend.trafo)
        self.trafo_choice = "rational"
        self.spacing = "equidistant"
        self.alpha = 0
        self.init_plot()

        with ui.row().classes('w-full no-wrap'):
            main_slider = ui.slider(min=0, max=self.slider_res, on_change=self.update_plot)
            main_slider.bind_value(self, "alpha")
            main_slider.tooltip("warp amount")

        with ((((ui.row())))):
            select_trafo = ui.select(list(self.available_trafos.keys()),
                                     value="rational",
                                     on_change=self.change_trafo)
            select_trafo.bind_value(self, "trafo_choice")
            select_trafo.tooltip("warp function")
            #select_trafo.props('popup-content-class="max-w-[200px]"')

            with ui.row().bind_visibility_from(self, "trafo_choice", backward=lambda x: x=="rational"):
                ui.label("pole/zero spacing")
                radio_spacing = ui.radio(["equidistant", "random"])
                radio_spacing.props('inline')
                radio_spacing.bind_value(self, "spacing")

        ui.run()

    def init_plot(self):
        self.main_plot = ui.matplotlib(figsize=(5, 5)).figure
        ax = self.main_plot.gca()
        self.pl1 = ax.plot(self.rend.hx, self.rend.hy, "b", dashes=[6, 2], linewidth=0.6)
        self.pl2 = ax.plot(self.rend.vx, self.rend.vy, "r", dashes=[2, 1], linewidth=0.6)
        ax.axis("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        plt.subplots_adjust(left=0.01,
                            right=0.99,
                            top=0.99,
                            bottom=0.01,
                            hspace = 0,
                            wspace = 0)
        plt.margins(0, 0)
        plt.tight_layout()

    def update_plot(self):
        self.rend.update(self.alpha/self.slider_res)
        with self.main_plot:
            for k, line in enumerate(self.pl1):
                line.set_data(self.rend.hx[k],
                              self.rend.hy[k])

            for k, line in enumerate(self.pl2):
                line.set_data(self.rend.vx[k],
                              self.rend.vy[k])

    def change_trafo(self):
        self.rend.trafo = self.available_trafos[self.trafo_choice]
        self.update_plot()


rend = renderer()
rendGUI = rendererGUI(rend)