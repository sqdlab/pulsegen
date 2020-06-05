'''
Created on May 25, 2013

@author: Matthias Baur
'''
import wx
import numpy as np
import copy

import matplotlib
matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from wxwidgets import Param, SliderGroup
        
# class PTPlotFrame(wx.Frame):
#     min_size = (600, 400)
# 
#     def __init__(self, parent=None, id_=-1, size=(900, 600),**kwargs):
#         self.size = self._size_constraint(size)
#         
#         super(PTPlotFrame, self).__init__(parent, id_, size=self.size, **kwargs)
#         
#         font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
#         font.SetPointSize(9)
#         
#         self._pt_plot_window = PTPlotWindow(self)
#         self._tstop_slider_group = SliderGroup(self, label="Stop Time (ns)",
#                         param=self._pt_plot_window.tstop_param)
#         self._twidth_slider_group = SliderGroup(self, label="Width (ns)",
#                         param=self._pt_plot_window.twidth_param)
#         self._pattern_slider_group = SliderGroup(self, label="Pattern Number",
#                         param=self._pt_plot_window.pattern_param)
#         
# #         self._toolbar = NavigationToolbar(self._pt_plot_window.canvas)
# 
#         grid_sizer = wx.GridBagSizer(3, 1)
#         grid_sizer.Add(self._tstop_slider_group.sizer, pos=(0, 0), flag=wx.EXPAND)
#         grid_sizer.Add(self._twidth_slider_group.sizer, pos=(1, 0), flag=wx.EXPAND)
#         grid_sizer.Add(self._pattern_slider_group.sizer, pos=(2, 0), flag=wx.EXPAND)
#         grid_sizer.AddGrowableCol(0)
#         
#         sizer = wx.BoxSizer(wx.VERTICAL)
# #         sizer.Add(self._toolbar, proportion=0, flag=wx.EXPAND)
#         sizer.Add(self._pt_plot_window, proportion=1, flag=wx.EXPAND | wx.GROW)
#         sizer.Add(grid_sizer, flag=wx.EXPAND | wx.LEFT | wx.RIGHT 
#                  | wx.BOTTOM | wx.TOP, border=10)
#         
#         self.SetSizer(sizer)
#         
#     def _size_constraint(self, size):
#         new_size = list(size)
#         if size[0] < size[1]:
#             return self.min_size
#         for i, value in enumerate(size):
#             if value < self.min_size[i]:
#                 new_size[i] = self.min_size[i]
#         return tuple(new_size)
#         
# class PTPlotWindow(wx.Window):
#     colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']
#     
#     def __init__(self, parent=None, id_=wx.ID_ANY, **kwargs):
#         super(PTPlotWindow, self).__init__(parent, id_, **kwargs)
#         
#         self.figure = Figure()
#         self.figure.patch.set_color('w')
#         self.canvas = FigureCanvas(self, -1, self.figure)
#         
#         self.tstop_param = Param(4000, 0, 6000)
#         self.twidth_param = Param(100, 0, 1000)
#         self.pattern_param = Param(0, 0, 1)
#         
#         self.tstop_param.attach(self)
#         self.twidth_param.attach(self)
#         self.pattern_param.attach(self)
#         
#         self._init_plots()
#         
#     def _init_plots(self):
#         self.waveform_plot = self.figure.add_subplot(211, 
#             xlabel=r"time (ns)", ylabel="amplitude")
#         self.marker_plot = self.figure.add_subplot(212, ylabel="amplitude")
#         
#         self.waveform_plot.set_ylim((-1.05, 1.05))
#         self.waveform_plot.set_xlim((0, 1))
#         self.waveform_plot.axhline(y=0, xmin=0, xmax=1, color='black')
#         
#         self.marker_plot.set_xticklabels([])
#         self.marker_plot.set_ylim((-0.1, 1.05))
#         self.marker_plot.set_xlim((0, 1))
#         self.Bind(wx.EVT_SIZE, self.size_handler)
#         
#     def size_handler(self, *args, **kwargs):
#         self.canvas.SetSize(self.GetSize())
#     
#     def draw(self):
#         pass
# #         x1, y1, x2, y2 = self.load_data(self)
#         
#     def draw_demo(self):
#         t = np.arange(0, 6000, 1)
#         y = np.sin(100*2 * np.pi * t/1000)
#         self.waveform_draw.plot(t, y)
#         self.waveform_plot.set_xlim((0,6000))
#         self.figure.tight_layout()
#         self.canvas.draw()
#         
#     def repaint(self):
#         self.figure.tight_layout()
#         self.canvas.draw()
#         
# class App(wx.App):
#     def OnInit(self):
#         self.frame1 = PTPlotFrame(title="Pulse Train Plot", size=(900, 600))
#         self.frame1.Show()
#         return True

def plot(seq, channels=None, markers=None, pattern=0):
    if channels is None:
        channels = list(range(len(seq.channels)))
    if markers is None:
        if(len(channels) 
            and len(seq.sequences)>channels[0]
            and len(seq.sequences[channels[0]])
        ):
            markers = list(range(len(seq.sequences[channels[0]].markers[0])))
        else:
            markers = []
    app = App(False)
    app.set_sequence(seq)
    app.plot(channels, markers, pattern)
    app.MainLoop()


class MainFrame(wx.Frame):
    '''
    classdocs
    
    '''
    # Make sequence attributes public
    title = "pulse train plot"
    tstop_default = 4000 + 10
    twidth_default = 100
    colors = ['b', 'r', 'g', 'c', 'm', 'y', 'k']

    def __init__(self, parent=None, id_=wx.ID_ANY, **kwargs):
        '''
        Constructor
        '''
        super(MainFrame, self).__init__(
            parent, id_, self.title, size=(900, 600), **kwargs)
        self.SetMinSize((600, 400))
        
        self._sequence = None
        self._channel_lines = {}
        self._marker_lines = {}
        self._tstop = self.tstop_default
        self._twidth = self.twidth_default
        self._tstart = self._tstop - self._twidth
        
        self._init_ui()
        
    def _init_ui(self):
        
        self._panel = wx.Panel(self)
        
        font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)
        font.SetPointSize(9)
        
        # setup the plot canvas
        self._figure = Figure((3, 2), dpi=80)
        self._figure.patch.set_color('w')
        self._plot_channels = self._figure.add_subplot(
                212, xlabel=r'time ($\mu s$)', ylabel="amplitude"
                )
        self._plot_markers = self._figure.add_subplot(211, ylabel="amplitude")
        self._plot_markers.set_xticklabels([])
#         self._figure.subplots_adjust()
        
        self._plot_channels.set_ylim((-1.05, 1.05))
        self._plot_channels.set_xlim((self._tstart, self._tstop))
        self._plot_channels.axhline(y=0, xmin=0, xmax=1, color='black')
        
        self._plot_markers.set_ylim((-0.1, 1.05))
        self._plot_markers.set_xlim((self._tstart, self._tstop))
        
        #self._plot_markers = self._figure.add_subplot(212)
        self._canvas = FigureCanvas(self._panel, -1, self._figure)
        self._toolbar = NavigationToolbar(self._canvas)
        
        # setup the controls
        self._text_ch = wx.StaticText(self._panel, label="Channels")
        self._ctrl_ch = wx.TextCtrl(self._panel, size=(100, -1), 
                                    style=wx.TE_PROCESS_ENTER)
        self._ctrl_ch.Bind(wx.EVT_TEXT_ENTER, self._ch_update)
        
        self._text_marker = wx.StaticText(self._panel, label="Markers")
        self._ctrl_marker = wx.TextCtrl(self._panel, size=(100, -1), 
                                        style=wx.TE_PROCESS_ENTER)
        self._ctrl_marker.Bind(wx.EVT_TEXT_ENTER, self._ch_update)
        
        self._text_stop = wx.StaticText(self._panel, label="Stop Time (ns)")
        
        self._slider_stop = wx.Slider(self._panel, value=self.tstop_default, 
                                      minValue=0, maxValue=6000)
#        self._slider_stop.SetMaxSize((400,-1))
        self._slider_stop.Bind(wx.EVT_SCROLL, self.on_slider_stop_scroll)
        
        self._text_width = wx.StaticText(self._panel, label="Width (ns)")
        
        self._slider_width = wx.Slider(self._panel, value=self.twidth_default, 
                                       minValue=1, maxValue=1000)
#         self._slider_width.SetMaxSize((400,-1))
        self._slider_width.Bind(wx.EVT_SCROLL, self.on_slider_width_scroll)
        
        self._button_stop_reset = wx.Button(self._panel, label="reset", size=(50, 30))
        self._button_stop_reset.Bind(wx.EVT_BUTTON, self._on_button_reset_stop)
        
        self._ctrl_width = wx.TextCtrl(self._panel, size=(60, -1), 
                                       value=str(self.twidth_default), 
                                       style=wx.TE_PROCESS_ENTER)
        self._ctrl_width.Bind(wx.EVT_TEXT_ENTER, self._on_ctrl_width)
        
        self._ctrl_stop = wx.TextCtrl(self._panel, size=(60, -1), 
                                      value=str(self.tstop_default), 
                                      style=wx.TE_PROCESS_ENTER)
        self._ctrl_stop.Bind(wx.EVT_TEXT_ENTER, self._on_ctrl_stop)
        
        self._button_width_reset = wx.Button(self._panel, label="reset", size=(50, 30))
        self._button_width_reset.Bind(wx.EVT_BUTTON, self._on_button_reset_width)
        
        self._text_pattern = wx.StaticText(self._panel, label="Pattern")
        
        self._slider_pattern = wx.Slider(self._panel, value=0, minValue=0, maxValue=1)
        self._slider_pattern.Bind(wx.EVT_SCROLL, self._on_slider_pattern)
        
        self._ctrl_pattern = wx.TextCtrl(self._panel, size=(60, -1), 
                                         value=str(0), style=wx.TE_PROCESS_ENTER)
        self._ctrl_pattern.Bind(wx.EVT_TEXT_ENTER, self._on_ctrl_pattern)
        
        # First row
        self._grid_sizer = wx.GridBagSizer(3, 6)
        
        self._grid_sizer.Add(self._text_ch, pos=(0, 0), flag=wx.ALIGN_RIGHT 
                       | wx.ALIGN_CENTER_VERTICAL)
        self._grid_sizer.Add(self._ctrl_ch, pos=(0, 1))
        self._grid_sizer.Add(self._text_stop, pos=(0, 2), flag=wx.LEFT | wx.ALIGN_RIGHT 
                       | wx.ALIGN_CENTER_VERTICAL, border=20)
        self._grid_sizer.Add(self._slider_stop, pos=(0, 3), flag=wx.EXPAND
                       | wx.ALIGN_RIGHT)
        self._grid_sizer.Add(self._ctrl_stop, pos=(0, 4), 
                             flag=wx.ALIGN_LEFT)
        self._grid_sizer.Add(self._button_stop_reset, pos=(0, 5), 
                             flag=wx.ALIGN_LEFT | wx.LEFT, border=10)

        # Second row
        self._grid_sizer.Add(self._text_marker, pos=(1, 0), flag=wx.ALIGN_RIGHT 
                       | wx.ALIGN_CENTER_VERTICAL)
        self._grid_sizer.Add(self._ctrl_marker, pos=(1, 1))
        self._grid_sizer.Add(self._text_width, pos=(1, 2), flag=wx.LEFT |wx.ALIGN_RIGHT 
                       | wx.ALIGN_CENTER_VERTICAL, border=20)
        self._grid_sizer.Add(self._slider_width, pos=(1, 3), flag=wx.EXPAND 
                       | wx.ALIGN_RIGHT)
        self._grid_sizer.Add(self._ctrl_width, pos=(1, 4), 
                             flag=wx.LEFT)
        self._grid_sizer.Add(self._button_width_reset, pos=(1, 5), 
                             flag=wx.ALIGN_LEFT | wx.LEFT, border=10)
        
        self._grid_sizer.Add(self._text_pattern, pos=(2, 2), flag=wx.LEFT 
                             | wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL, border=20)
        self._grid_sizer.Add(self._slider_pattern, pos=(2, 3), flag=wx.EXPAND)
        self._grid_sizer.Add(self._ctrl_pattern, pos=(2, 4), flag=wx.LEFT)
        
        self._grid_sizer.AddGrowableCol(3, 1)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self._toolbar, proportion=0, flag=wx.LEFT | wx.EXPAND)
        vbox.Add(self._canvas, proportion=1, flag=wx.EXPAND)
        vbox.Add(self._grid_sizer, flag=wx.EXPAND | wx.LEFT | wx.RIGHT 
                 | wx.BOTTOM | wx.TOP, border=10)
        self._panel.SetSizer(vbox)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self._panel, proportion=1, flag=wx.EXPAND)
        
        self.SetSizer(sizer)
        
        self.Centre()
        self.Show()
        
    def _create_channel_line(self, ch_idx, pat_idx):
        key = "{0}, {1}".format(ch_idx, pat_idx)
        if key in self._channel_lines:
            return self._channel_lines[key]
        
        channel = self._many_awg_seq.channels[ch_idx]
        tmax = channel.pattern_length * 1.0e9
        Ts = 1/channel.sampling_freq * 1.0e9
        
        times = np.arange(0, tmax, Ts)
        sequence = self._many_awg_seq.sampled_sequences[ch_idx]
        ypoints = sequence.segments[pat_idx].waveform
        
        self._channel_lines[key] = Line2D(times, ypoints, marker='o', 
                                              markersize=3)
        return self._channel_lines[key]
    
    def _create_marker_line(self, ch_idx, m_idx, pat_idx):
        key = "{0}, {1}, {2}".format(ch_idx, m_idx, pat_idx)
        if key in self._marker_lines:
            return self._marker_lines[key]
        
        channel = self._many_awg_seq.channels[ch_idx]
        tmax = channel.pattern_length * 1.0e9
        sampling_freq = channel.sampling_freq/channel.awg.granularity
        Ts = 1.0e9/sampling_freq
        
        sequence = self._many_awg_seq.sampled_sequences[ch_idx]
        ypoints = sequence.segments[pat_idx].markers[m_idx]
        times = Ts*np.arange(0, len(ypoints))
        
        self._marker_lines[key] = Line2D(times, ypoints)
        return self._marker_lines[key]
    
    def plot(self, channels, markers, pattern):
        
        # update the controls if the function is called by the user
        self._ctrl_ch.SetValue(','.join(map(str, channels)))
        self._ctrl_marker.SetValue(','.join(map(str, markers)))
        self._ctrl_pattern.SetValue(str(pattern))
        self._slider_pattern.SetValue(pattern)
        
        if self._many_awg_seq is None:
            return
        
        self._plot_channels.lines = []
        self._plot_markers.lines = []
        
        if not channels:
            return
        
        number_of_markers = len(channels) * len(markers)
        if number_of_markers:
            marker_amp = 1.0/number_of_markers
         
        marker_idx = 0
        for ch_idx, channel in enumerate(channels):
            current_line = self._create_channel_line(channel, pattern)
            current_line.set_color(self.colors[ch_idx])
            self._plot_channels.add_line(current_line)
            
            for marker in markers:
                current_line = copy.copy(
                    self._create_marker_line(channel, marker, pattern))
                x, y = current_line.get_data()
                y = y*0.9*marker_amp + marker_idx*marker_amp
                current_line.set_data(x, y)
                current_line.set_color(self.colors[ch_idx])
                self._plot_markers.add_line(current_line)
                
                marker_idx += 1
        
        self._figure.tight_layout()
        self._figure.canvas.draw_idle()
        
    def plot_demo(self):
        t = np.arange(0, 6000, 1)
        y = np.sin(100*2 * np.pi * t/1000)
        self._plot_channels.plot(t, y)
        self._figure.tight_layout()
        
    def _parse_ctrl_str(self, ctrl_object):
        ch_str = ctrl_object.GetValue()
        if ch_str:
            return [int(ch) for ch in ch_str.split(',')]
        else:
            return []
        
    def _ch_update(self, event):
        ch_list = self._parse_ctrl_str(self._ctrl_ch)
        marker_list = self._parse_ctrl_str(self._ctrl_marker)
        pattern = self._slider_pattern.GetValue()
        self.plot(ch_list, marker_list, pattern)
    
    def set_sequence(self, sequence_obj):
        self._many_awg_seq = sequence_obj
        sequences = self._many_awg_seq.sequences
        channels = self._many_awg_seq.channels
        self._slider_pattern.SetMax(max([len(seq) for seq in sequences]) - 1)
        
        self._tstop = int(channels[0].fixed_points[0]*1.0e9 + 10)
        self._new_tstop_twidth()
        
        pattern_length = int(channels[0].pattern_length * 1.0e9)
        self._slider_stop.SetMax(pattern_length)
    
    def _new_tstop_twidth(self):
        self._tstart = self._tstop - self._twidth
        self._ctrl_stop.SetValue(str(self._tstop))
        self._ctrl_width.SetValue(str(self._twidth))
        self._plot_channels.set_xlim((self._tstart, self._tstop))
        self._plot_markers.set_xlim((self._tstart, self._tstop))
        self._figure.canvas.draw_idle()
    
    def on_slider_stop_scroll(self, event):
        self._tstop = event.GetPosition()
        self._ctrl_stop.SetValue(str(self._tstop))
        self._new_tstop_twidth()
    
    def on_slider_width_scroll(self, event):
        self._twidth = event.GetPosition()
        self._ctrl_width.SetValue(str(self._twidth))
        self._new_tstop_twidth()
        
    def _on_slider_pattern(self, event):
        self._ch_update(event)
        
    def _on_button_reset_stop(self, event):
        self._tstop = self.tstop_default
        self._ctrl_stop.SetValue(str(self._tstop))
        self._slider_stop.SetValue(self._tstop)
        self._new_tstop_twidth()
        
    def _on_button_reset_width(self, event):
        self._twidth = self.twidth_default
        self._ctrl_width.SetValue(str(self._twidth))
        self._slider_width.SetValue(self._twidth)
        self._new_tstop_twidth()
        
    def _on_ctrl_width(self, event):
        value = int(self._ctrl_width.GetValue())
        self._twidth = value
        self._slider_width.SetValue(value)
        self._new_tstop_twidth()
    
    def _on_ctrl_stop(self, event):
        value = int(self._ctrl_stop.GetValue())
        self._tstop = value
        self._slider_stop.SetValue(value)
        self._new_tstop_twidth()
        
    def _on_ctrl_pattern(self, event):
        value = int(self._ctrl_pattern.GetValue())
        self._slider_pattern.SetValue(value)
        self._ch_update(event)


class App(wx.App):
    def OnInit(self):
        self.frame = MainFrame()
        self.frame.Show()
        return True
    
    def plot(self, channels, markers, pattern):
        self.frame.plot(channels, markers, pattern)
    
    def set_sequence(self, sequence_obj):
        self.frame.set_sequence(sequence_obj)


def main():
    app = App(False)
    app.frame.plot_demo()
    app.MainLoop()

if __name__ == '__main__':
    main()