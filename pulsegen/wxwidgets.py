import wx.lib.newevent
import cPickle

DragListEvent, EVT_DRAGLIST = wx.lib.newevent.NewEvent()

class DragListCtrl(wx.ListCtrl):
    '''
    Taken and modified from http://wiki.wxpython.org/ListControls
    '''
    def __init__(self, parent, *args, **kwargs):
        super(DragListCtrl, self).__init__(parent, *args, **kwargs)

        self.Bind(wx.EVT_LIST_BEGIN_DRAG, self._startDrag)

        dt = ListDrop(self)
        self.SetDropTarget(dt)

    def getItemInfo(self, idx):
        """Collect all relevant data of a listitem, and put it in a list"""
        l = []
        l.append(idx) # We need the original index, so it is easier to eventualy delete it
        l.append(self.GetItemData(idx)) # Itemdata
        l.append(self.GetItemText(idx)) # Text first column
        for i in range(1, self.GetColumnCount()): # Possible extra columns
            l.append(self.GetItem(idx, i).GetText())
        return l

    def _startDrag(self, e):
        """ Put together a data object for drag-and-drop _from_ this list. """
        l = []
        idx = -1
        while True: # find all the selected items and put them in a list
            idx = self.GetNextItem(idx, wx.LIST_NEXT_ALL, wx.LIST_STATE_SELECTED)
            if idx == -1:
                break
            l.append(self.getItemInfo(idx))

        # Pickle the items list.
        itemdata = cPickle.dumps(l, 1)
        # create our own data format and use it in a
        # custom data object
        ldata = wx.CustomDataObject("ListCtrlItems")
        ldata.SetData(itemdata)
        # Now make a data object for the  item list.
        data = wx.DataObjectComposite()
        data.Add(ldata)

        # Create drop source and begin drag-and-drop.
        dropSource = wx.DropSource(self)
        dropSource.SetData(data)
        res = dropSource.DoDragDrop(flags=wx.Drag_DefaultMove)

        # If move, we want to remove the item from this list.
        if res == wx.DragMove:
            # It's possible we are dragging/dropping from this list to this list.  In which case, the
            # index we are removing may have changed...

            # Find correct position.
            l.reverse() # Delete all the items, starting with the last item
            for i in l:
                pos = self.FindItem(i[0], i[2])
                self.DeleteItem(pos)

    def _insert(self, x, y, seq):
        """ Insert text at given x, y coordinates --- used with drag-and-drop. """

        # Find insertion point.
        index, flags = self.HitTest((x, y))

        if index == wx.NOT_FOUND: # not clicked on an item
            if flags & (wx.LIST_HITTEST_NOWHERE|wx.LIST_HITTEST_ABOVE|wx.LIST_HITTEST_BELOW): # empty list or below last item
                index = self.GetItemCount() # append to end of list
            elif self.GetItemCount() > 0:
                if y <= self.GetItemRect(0).y: # clicked just above first item
                    index = 0 # append to top of list
                else:
                    index = self.GetItemCount() + 1 # append to end of list
        else: # clicked on an item
            # Get bounding rectangle for the item the user is dropping over.
            rect = self.GetItemRect(index)

            # If the user is dropping into the lower half of the rect, we want to insert _after_ this item.
            # Correct for the fact that there may be a heading involved
            if y > rect.y - self.GetItemRect(0).y + rect.height/2:
                index += 1

        for i in seq: # insert the item data
            idx = self.InsertStringItem(index, i[2])
            self.SetItemData(idx, i[1])
            for j in range(1, self.GetColumnCount()):
                try: # Target list can have more columns than source
                    self.SetStringItem(idx, j, i[2+j])
                except:
                    pass # ignore the extra columns
            index += 1

class ListDrop(wx.PyDropTarget):
    """ Drop target for simple lists. """

    def __init__(self, source):
        """ Arguments:
         - source: source listctrl.
        """
        super(ListDrop, self).__init__()

        self.dv = source

        # specify the type of data we will accept
        self.data = wx.CustomDataObject("ListCtrlItems")
        self.SetDataObject(self.data)

    # Called when OnDrop returns True.  We need to get the data and
    # do something with it.
    def OnData(self, x, y, d):
        # copy the data from the drag source to our data object
        if self.GetData():
            # convert it back to a list and give it to the viewer
            ldata = self.data.GetData()
            l = cPickle.loads(ldata)
            self.dv._insert(x, y, l)
            evt = DragListEvent()
            wx.PostEvent(self.dv, evt)

        # what is returned signals the source what to do
        # with the original data (move, copy, etc.)  In this
        # case we just return the suggested value given to us.
        return d

class Param(object):
    """
    The idea of the "Param" class is that some parameter in the GUI may have
    several knobs that both control it and reflect the parameter's state, e.g.
    a slider, text, and dragging can all change the value of the frequency in
    the waveform of this example.
    The class allows a cleaner way to update/"feedback" to the other knobs when
    one is being changed.  Also, this class handles min/max constraints for all
    the knobs.
    Idea - knob list - in "set" method, knob object is passed as well
      - the other knobs in the knob list have a "set" method which gets
        called for the others.
    """
    def __init__(self, initialValue=None, minimum=0, maximum=1):
        self.minimum = minimum
        self.maximum = maximum
        if initialValue != self.constrain(initialValue):
            raise ValueError('illegal initial value')
        self.value = initialValue
        self.knobs = []

    def attach(self, knob):
        self.knobs += [knob]

    def set(self, value, knob=None):
        self.value = self.constrain(value)
        for feedbackKnob in self.knobs:
            if feedbackKnob != knob:
                feedbackKnob.set_knob(self.value)
        return self.value

    def constrain(self, value):
        if value <= self.minimum:
            value = self.minimum
        if value >= self.maximum:
            value = self.maximum
        return value

class SliderGroup(wx.Panel):
    # TODO: Fix flickering of slider
    def __init__(self, parent, label, param, id_=wx.ID_ANY):
        super(SliderGroup, self).__init__(parent, id=id_)
        self.id = id_
        self._label = wx.StaticText(self, label=label)
        self.text = wx.TextCtrl(self, -1, style=wx.TE_PROCESS_ENTER)
        self.slider = wx.Slider(self, -1)
        #self.slider.SetMax(param.maximum)
        self.set_knob(param.value)
        
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self._label, proportion=0, flag=wx.EXPAND 
                  | wx.ALIGN_CENTER | wx.ALL, border=5)
        sizer.Add(self.slider, proportion=1, flag=wx.EXPAND)
        sizer.Add(self.text, proportion=0, flag=wx.EXPAND 
                    | wx.ALIGN_CENTER | wx.ALL, border=5)
        self.sizer = sizer
        self.SetSizer(sizer)
        
        self.slider.Bind(wx.EVT_SLIDER, self._slider_handler)
        self.text.Bind(wx.EVT_TEXT_ENTER, self._text_handler)
        
        self._param = param
        self._param.attach(self)
        
    def _slider_handler(self, event):
        value = event.GetInt()
        self._param.set(value)
        event.Skip(True)
        
    def _text_handler(self, event):
        value = int(self.text.GetValue())
        self._param.set(value)
        event.Skip(True)
        
    def set_knob(self, value):
        self.text.SetValue(str(value))
        self.text.Refresh()
        self.text.Update()
        self.slider.SetValue(value)
        
    def disable(self):
        self.enable(False)
        
    def enable(self, enable=True):
        for ctrl in [self.text, self.slider, self._label]:
            ctrl.Enable(enable)
        
    def delete(self):
        self.slider.Unbind(wx.EVT_SLIDER)
        self.text.Unbind(wx.EVT_TEXT_ENTER)
        # TODO: replace with self.DestroyChildren()?
        self.text.Destroy()
        self.slider.Destroy()
        self._label.Destroy()
        del self