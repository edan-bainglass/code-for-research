"""
GUI neutral widgets
===================

Widgets that are designed to work for any of the GUI backends.
All of these widgets require you to predefine a :class:`matplotlib.axes.Axes`
instance and pass that as the first arg.  matplotlib doesn't try to
be too smart with respect to layout -- you will have to figure out how
wide and tall you want your Axes to be to accommodate your widget.
"""

import numpy as np
import re

from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt


class Widget(object):
    """
    Abstract base class for GUI neutral widgets
    """
    drawon = True
    eventson = True
    _active = True

    def set_active(self, active):
        """Set whether the widget is active.
        """
        self._active = active

    def get_active(self):
        """Get whether the widget is active.
        """
        return self._active

    # set_active is overridden by SelectorWidgets.
    active = property(get_active,
                      lambda self, active: self.set_active(active),
                      doc="Is the widget active?")

    def ignore(self, event):
        """Return True if event should be ignored.

        This method (or a version of it) should be called at the beginning
        of any event callback.
        """
        return not self.active


class AxesWidget(Widget):
    """Widget that is connected to a single
    :class:`~matplotlib.axes.Axes`.

    To guarantee that the widget remains responsive and not garbage-collected,
    a reference to the object should be maintained by the user.

    This is necessary because the callback registry
    maintains only weak-refs to the functions, which are member
    functions of the widget.  If there are no references to the widget
    object it may be garbage collected which will disconnect the
    callbacks.

    Attributes:

    *ax* : :class:`~matplotlib.axes.Axes`
        The parent axes for the widget
    *canvas* : :class:`~matplotlib.backend_bases.FigureCanvasBase` subclass
        The parent figure canvas for the widget.
    *active* : bool
        If False, the widget does not respond to events.
    """
    def __init__(self, ax):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.cids = []

    def connect_event(self, event, callback):
        """Connect callback with an event.

        This should be used in lieu of `figure.canvas.mpl_connect` since this
        function stores callback ids for later clean up.
        """
        cid = self.canvas.mpl_connect(event, callback)
        self.cids.append(cid)

    def disconnect_events(self):
        """Disconnect all events created by this widget."""
        for c in self.cids:
            self.canvas.mpl_disconnect(c)


class Slider(AxesWidget):
    """
    A slider representing a floating point range.

    Create a slider from *valmin* to *valmax* in axes *ax*. For the slider to
    remain responsive you must maintain a reference to it. Call
    :meth:`on_changed` to connect to the slider event.

    Attributes
    ----------
    val : float
        Slider value.
    """
    def __init__(self,
                 ax,
                 label,
                 valmin,
                 valmax,
                 valinit=0.5,
                 valfmt='%1.2f',
                 closedmin=True,
                 closedmax=True,
                 slidermin=None,
                 slidermax=None,
                 dragging=True,
                 valstep=None,
                 **kwargs):
        """
        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in.

        label : str
            Slider label.

        valmin : float
            The minimum value of the slider.

        valmax : float
            The maximum value of the slider.

        valinit : float, optional, default: 0.5
            The slider initial position.

        valfmt : str, optional, default: "%1.2f"
            Used to format the slider value, fprint format string.

        closedmin : bool, optional, default: True
            Indicate whether the slider interval is closed on the bottom.

        closedmax : bool, optional, default: True
            Indicate whether the slider interval is closed on the top.

        slidermin : Slider, optional, default: None
            Do not allow the current slider to have a value less than
            the value of the Slider `slidermin`.

        slidermax : Slider, optional, default: None
            Do not allow the current slider to have a value greater than
            the value of the Slider `slidermax`.

        dragging : bool, optional, default: True
            If True the slider can be dragged by the mouse.

        valstep : float, optional, default: None
            If given, the slider will snap to multiples of `valstep`.

        Notes
        -----
        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`~matplotlib.patches.Rectangle` that draws the slider
        knob.  See the :class:`~matplotlib.patches.Rectangle` documentation for
        valid property names (e.g., `facecolor`, `edgecolor`, `alpha`).
        """
        AxesWidget.__init__(self, ax)

        if slidermin is not None and not hasattr(slidermin, 'val'):
            raise ValueError("Argument slidermin ({}) has no 'val'".format(
                type(slidermin)))
        if slidermax is not None and not hasattr(slidermax, 'val'):
            raise ValueError("Argument slidermax ({}) has no 'val'".format(
                type(slidermax)))
        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False
        self.valmin = valmin
        self.valmax = valmax
        self.valstep = valstep
        valinit = self._value_in_bounds(valinit)
        if valinit is None:
            valinit = valmin
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axvspan(valmin, valinit, 0, 1, **kwargs)
        self.vline = ax.axvline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_yticks([])
        ax.set_xlim((valmin, valmax))
        ax.set_xticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(-0.02,
                             0.5,
                             label,
                             transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='right')

        self.valtext = ax.text(1.02,
                               0.5,
                               valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='left')

        self.cnt = 0
        self.observers = {}

        self.set_val(valinit)

    def _value_in_bounds(self, val):
        """ Makes sure self.val is with given bounds."""
        if self.valstep:
            val = np.round((val - self.valmin) / self.valstep) * self.valstep
            val += self.valmin

        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val
        return val

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif (
            (event.name == 'button_release_event') or
            (event.name == 'button_press_event' and event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return
        val = self._value_in_bounds(event.xdata)
        if val not in [None, self.val]:
            self.set_val(val)

    def set_val(self, val):
        """
        Set slider value to *val*

        Parameters
        ----------
        val : float
        """
        xy = self.poly.xy
        xy[2] = val, 1
        xy[3] = val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed call *func* with the new
        slider value

        Parameters
        ----------
        func : callable
            Function to call when slider is changed.
            The function must accept a single float as its arguments.

        Returns
        -------
        cid : int
            Connection id (which can be used to disconnect *func*)
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """
        Remove the observer with connection id *cid*

        Parameters
        ----------
        cid : int
            Connection id of the observer to be removed
        """
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """Reset the slider to the initial value"""
        if self.val != self.valinit:
            self.set_val(self.valinit)


class CheckButtons(AxesWidget):
    """
    A GUI neutral set of check buttons.

    For the check buttons to remain responsive you must keep a
    reference to this object.

    The following attributes are exposed

     *ax*
        The :class:`matplotlib.axes.Axes` instance the buttons are
        located in

     *labels*
        List of :class:`matplotlib.text.Text` instances

     *lines*
        List of (line1, line2) tuples for the x's in the check boxes.
        These lines exist for each box, but have ``set_visible(False)``
        when its box is not checked.

     *rectangles*
        List of :class:`matplotlib.patches.Rectangle` instances

    Connect to the CheckButtons with the :meth:`on_clicked` method
    """
    def __init__(self, ax, labels, actives=None):
        """
        Add check buttons to :class:`matplotlib.axes.Axes` instance *ax*

        Parameters
        ----------
        ax : `~matplotlib.axes.Axes`
            The parent axes for the widget.

        labels : List[str]
            The labels of the check buttons.

        actives : List[bool], optional
            The initial check states of the buttons. The list must have the
            same length as *labels*. If not given, all buttons are unchecked.
        """
        AxesWidget.__init__(self, ax)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_navigate(False)

        if actives is None:
            actives = [False] * len(labels)

        if len(labels) > 1:
            dy = 1. / (len(labels) + 1)
            ys = np.linspace(1 - dy, dy, len(labels))
        else:
            dy = 0.25
            ys = [0.5]

        axcolor = ax.get_facecolor()

        self.labels = []
        self.lines = []
        self.rectangles = []

        lineparams = {
            'color': 'k',
            'linewidth': 1.25,
            'transform': ax.transAxes,
            'solid_capstyle': 'butt'
        }
        for y, label, active in zip(ys, labels, actives):
            t = ax.text(0.25,
                        y,
                        label,
                        transform=ax.transAxes,
                        horizontalalignment='left',
                        verticalalignment='center')

            w, h = dy / 2, dy / 2
            x, y = 0.05, y - h / 2

            p = Rectangle(xy=(x, y),
                          width=w,
                          height=h,
                          edgecolor='black',
                          facecolor=axcolor,
                          transform=ax.transAxes)

            l1 = Line2D([x, x + w], [y + h, y], **lineparams)
            l2 = Line2D([x, x + w], [y, y + h], **lineparams)

            l1.set_visible(active)
            l2.set_visible(active)
            self.labels.append(t)
            self.rectangles.append(p)
            self.lines.append((l1, l2))
            ax.add_patch(p)
            ax.add_line(l1)
            ax.add_line(l2)

        self.connect_event('button_press_event', self._clicked)

        self.cnt = 0
        self.observers = {}

    def _clicked(self, event):
        if self.ignore(event) or event.button != 1 or event.inaxes != self.ax:
            return
        for i, (p, t) in enumerate(zip(self.rectangles, self.labels)):
            if (t.get_window_extent().contains(event.x, event.y)
                    or p.get_window_extent().contains(event.x, event.y)):
                self.set_active(i)
                break

    def set_active(self, index):
        """
        Directly (de)activate a check button by index.

        *index* is an index into the original label list
            that this object was constructed with.
            Raises ValueError if *index* is invalid.

        Callbacks will be triggered if :attr:`eventson` is True.

        """
        if 0 > index >= len(self.labels):
            raise ValueError("Invalid CheckButton index: %d" % index)

        l1, l2 = self.lines[index]
        l1.set_visible(not l1.get_visible())
        l2.set_visible(not l2.get_visible())

        if self.drawon:
            self.ax.figure.canvas.draw()

        if not self.eventson:
            return
        for cid, func in self.observers.items():
            func(self.labels[index].get_text())

    def get_status(self):
        """
        returns a tuple of the status (True/False) of all of the check buttons
        """
        return [l1.get_visible() for (l1, l2) in self.lines]

    def on_clicked(self, func):
        """
        When the button is clicked, call *func* with button label

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass


class PremiumCheckButtons(CheckButtons, AxesWidget):
    def __init__(self, ax, lines, actives, linecolor="k", showedge=True, **kw):
        AxesWidget.__init__(self, ax)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_navigate(False)
        if not showedge:
            ax.axis("off")
        linekw = {'solid_capstyle': 'butt', "color": linecolor}
        colors = [l.get_color() for l in lines]
        i = []

        class Handler(object):
            def legend_artist(self, legend, orig_handle, fontsize, handlebox):
                x0, y0 = handlebox.xdescent, handlebox.ydescent
                height = handlebox.height
                self.line1 = plt.Line2D([x0, x0 + height], [y0, y0 + height],
                                        **linekw)
                self.line2 = plt.Line2D([x0, x0 + height], [y0 + height, y0],
                                        **linekw)
                self.rect = plt.Rectangle((x0, y0),
                                          height,
                                          height,
                                          edgecolor="k",
                                          fill=False)
                self.line = plt.Line2D([x0 + 2 * height, x0 + 3.2 * height],
                                       [y0 + 0.5 * height, y0 + 0.5 * height],
                                       color=colors[len(i)],
                                       alpha=0.75,
                                       linewidth=2)
                handlebox.add_artist(self.rect)
                handlebox.add_artist(self.line)
                handlebox.add_artist(self.line1)
                handlebox.add_artist(self.line2)
                i.append(1)
                return [self.line1, self.line2, self.rect]

        labels = [re.sub('([0-9]+)', r'$_{\1}$', l.get_label()) for l in lines]
        self.box = ax.legend(handles=[object() for l in lines],
                             labels=labels,
                             handler_map={object: Handler()},
                             **kw)

        self.lines = [(h[0], h[1]) for h in self.box.legendHandles]
        self.rectangles = [h[2] for h in self.box.legendHandles]
        self.labels = self.box.texts

        for i, (l1, l2) in enumerate(self.lines):
            l1.set_visible(actives[i])
            l2.set_visible(actives[i])

        self.connect_event('button_press_event', self._clicked)

        self.cnt = 0
        self.observers = {}