#!/usr/bin/python3
# -*- coding: utf-8 -*-
# -*- coding: iso-8859-15 -*-
# -*- coding: latin-1 -*-
# -*- coding: ascii -*-
from __future__ import print_function

import time, sys
from threading import *
import wx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

from mpl_toolkits.mplot3d import Axes3D

from xrr import XRR, _XLayer
from xrr_wxg import XRR_Frame



class event_socket(wx.PyEvent):
    def __init__(self, wxgui=None, key=1, name='', func=None):
        wx.PyEvent.__init__(self)  
        self.wxgui = wxgui
        self.key = key
        self.name = name
        self.func = func
        self.EVTID = wx.NewId()
        self.SetEventType(self.EVTID)
        self.wxgui.Connect(-1, -1, self.EVTID, self.func )
    def write(self, text): 
        self.data = text
        wx.PostEvent(self.wxgui, self)    
    def flush(self):       pass;
    
    
class event_socket_dict():
    def __init__(self, wxgui=None):
        self.wxgui = wxgui
        self.dict = {} 
        
    def add(self, key=1, name='', func=None):
        self.dict.update( {key : event_socket(self.wxgui, key, name, func) } )
        
    def post(self, key=1, data={}):
        self.dict[key].write(data)

        
class XRR_wx(XRR):
    # ~ def __init__(self, *args, **kwds):
    # ~     super().__init__()
    
    def __init__(self, event_dict):
        print('XRR_wx: __init__()')
        self.event_dict = event_dict
        XRR.__init__(self)
        
    def post_event(self, **kw):
    #=================================
        # ~ print(key, data, message, end='', flush=True)
        # ~ # can be redifined when the class is overloaded
        # ~ self.pyqtsig_val.emit(data)
        m = {0:'L-M', 1:'D-E'}
        d = {'@':'XRR', 'fit':m[self.fitmode], 'fitn':self.fitn}
        d.update(kw)
        # ~ print(d)
        key = kw['key']
        self.event_dict.post( key, data=d)
        return
    

        
# Thread class that executes processing
class WorkerThread(Thread):
    """Worker Thread Class."""

    def __init__(self, event_dict):
        """Init Worker Thread Class."""
        print('>>#worker>: __init__()')
        Thread.__init__(self)
        self.event_dict = event_dict
        self.wxgui = self.event_dict.wxgui
        self.np = 0
        # ~ self.start() # better from the wxgui
        
    def post_event_wk(self, key=1, msg=''):
        self.np += 1
        d = {'@':'worker', 'key':key, 'msg':msg}
        print(d, flush= True)
        self.event_dict.post(key, data=d)

    def run(self):
        self.post_event_wk(key=1, msg='run start' )
        self.wxgui.xrr.layer_print()
        self.wxgui.xrr.fit()
        self.post_event_wk(key=0, msg='run end' )

    def want_abort(self):
        """abort worker thread. request"""
        # Method for use by wxgui to request abort
        self.wxgui.xrr.abort = 1  # here ... just pass it to the xrr
        self.post_event_wk(key=-1, msg='abort' )



        
class Plot_Notebook(wx.Panel):
    def __init__(self, parent, id=-1):
        wx.Panel.__init__(self, parent, id=id)
        self.notebook = wx.Notebook(self, wx.ID_ANY)
        sizer = wx.BoxSizer()
        sizer.Add(self.notebook, 1, wx.EXPAND)
        self.SetSizer(sizer)
        self.pagelist = []

    def add(self, name="plot"):
        page = Plot_Panel(self.notebook)
        self.notebook.AddPage(page, name)
        self.pagelist.append(page)
        return page
        
        
class Plot_Panel(wx.Panel):
    def __init__(self, parent, id=-1, dpi=None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)  
        self.figure = mpl.figure.Figure(dpi=dpi)    #, figsize=(2, 2)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.toolbar = NavigationToolbar(self.canvas)
        self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        # ~ sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)
        
class Plot_Panel0(wx.Panel):
    def __init__(self, parent, id=-1, dpi=None, **kwargs):
        wx.Panel.__init__(self, parent, id=id, **kwargs)  
        self.figure = mpl.figure.Figure(dpi=dpi)    #, figsize=(2, 2)
        self.canvas = FigureCanvas(self, -1, self.figure)
        # ~ self.toolbar = NavigationToolbar(self.canvas)
        # ~ self.toolbar.Realize()

        sizer = wx.BoxSizer(wx.VERTICAL)
        # ~ sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        # ~ sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)


# GUI Frame class that spins off the worker thread
# ~ class MainFrame(wx.Frame):
class MainFrame(XRR_Frame):
    """Class MainFrame."""
        
    # ~ def __init__(self):  #, *args, **kwds):
        # ~ super().__init__()
    def __init__(self, *args, **kwds):
        XRR_Frame.__init__(self, *args, **kwds)
        
        self.SetTitle("X-ray reflectivity - Parrat algorithm")
        
        self.contentNotSaved = True
        
        self.tx = self.text_ctrl_log        # ~ wx.TextCtrl(self, wx.NewId(), style=wx.TE_RICH)

        self.Bind(wx.EVT_BUTTON, self.OnStart, self.button_fit)
        
        self.EvtD_xrrfit = event_socket_dict(wxgui=self)
        self.EvtD_xrrfit.add(key=1,  name='xrr_run', func=self.OnEvt_xrr_run)
        self.EvtD_xrrfit.add(key=0,  name='xrr_end', func=self.OnEvt_xrr_end)
        self.EvtD_xrrfit.add(key=-1, name='xrr_abt', func=self.OnEvt_xrr_abt)
        self.xrr = XRR_wx(self.EvtD_xrrfit)
        
        self.EvtD_worker = event_socket_dict(wxgui=self)
        self.EvtD_worker.add(key=1,  name='wrk_run', func=self.OnEvt_wrk_all)
        self.EvtD_worker.add(key=0,  name='wrk_end', func=self.OnEvt_wrk_all)
        self.EvtD_worker.add(key=-1, name='wrk_abt', func=self.OnEvt_wrk_all)
        # ~ self.worker = WorkerThread(self.EvtD_worker)
        
        '''
        # Set up event handler for any worker thread results
        # And indicate we don't have a worker thread yet
        '''
        self.worker = None
        self.nrun =0
        
        self.layer_to_grid()
        
        self.plotter1 = Plot_Notebook(self.panel_1)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.plotter1, -1, wx.ALL|wx.EXPAND, 5)
        self.panel_1.SetSizer(sizer)
        self.plotter1.add('Int vs 2tht')
        self.plotter1.add('Int*Q^4 vs Q')
        self.plotter1.add('Layer profile')
        self.data_plot1()
        self.data_plot2()
        
        self.plot3 = Plot_Panel0(self.panel_2)
        self.plot3.figure.add_subplot(111)
        sizer3_ = wx.BoxSizer(wx.VERTICAL)
        sizer3_.Add(self.plot3, -1, wx.ALL|wx.EXPAND, 5)
        self.panel_2.SetSizer(sizer3_)
        
        # ~ self.xrr.layer_profile()
        self.layer_plot()
        
        sys.stdout = event_socket(wxgui=self, key=0, name='stdout', func=self.OnPrint)
        sys.stderr = event_socket(wxgui=self, key=1, name='stderr', func=self.OnPrint)
        
        print('@wxqui>: init')
        print('@wxqui>: init - finished')
        
    def __del__(self):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        
    def OnPrint(self, event):
        name = event.name
        if name is 'stdout':
            font = wx.Font( wx.FontInfo(8).Family(wx.TELETYPE) )
            self.tx.SetForegroundColour(wx.BLACK)
        if name is 'stderr':
            font = wx.Font( wx.FontInfo(8).Bold().Family(wx.TELETYPE) )
            self.tx.SetForegroundColour(wx.RED)
            
        self.text_ctrl_log.SetFont(font)
        self.text_ctrl_log.write(event.data)
        
        # ~ self.text_ctrl_log.write(str(event.data))
        # ~ self.text_ctrl_log.WriteText(str(event.data))
        # ~ wx.CallAfter(self.text_ctrl_log.write, str(event.data))    #thread safe
        # ~ wx.CallAfter(self.text_ctrl_log.WriteText, str(event.data))    #thread safe
        # ~ wx.Yield()
        

    def OnStart(self, event):
        # Trigger the worker thread unless it's already busy
        if self.worker == None:
            """Start Computation."""
            self.nrun =0
            self.tx.Clear()
            print('@wxqui> Creating worker and start')
            self.worker = WorkerThread(self.EvtD_worker)

            self.button_fit.SetLabel(' Abort ')
            self.worker.start()
        else:
            if self.worker.is_alive():
                """Abort Computation."""
                print('@wxqui>: Send worker.want_abort()')
                self.worker.want_abort()
                self.button_fit.SetLabel(' .... ')
            else:  
                """Reseting Computation."""
                self.evget_w_end
            
    def OnEvt_wrk_all(self, event):
        print( '@', event.name, '<<', event.data, )
        
    def OnEvt_xrr_update(self):
        self.layer_to_grid()
        # ~ self.xrr.layer_profile()
        # ~ self.layer_plot()
        self.data_plot1()
        self.data_plot2()
        self.data_plot3()
            
    def OnEvt_xrr_end(self, event):
        print( '@', event.name, '<<', event.data, )
        self.worker = None
        self.button_fit.SetLabel('Fit Data')
        self.OnEvt_xrr_update()
        self.xrr.layer_profile()
        self.layer_plot()
        
    def OnEvt_xrr_abt(self, event):
        print( '@', event.name, '<<', event.data, )
        
    def OnEvt_xrr_run(self, event):
        print( '@', event.name, '<<', event.data, )
        mx = self.slider_fit.GetMax()
        nf = event.data['fitn']
        if (nf  % 100) == 0:
            print('==========================')
            self.OnEvt_xrr_update()
        if  nf > mx:
            self.slider_fit.SetMax(mx+1000)
        self.slider_fit.SetValue(self.xrr.fitn)


        
    ################################################################
        
    # ~ def plot_profile(self):
    def layer_plot(self):
        pg = self.plotter1.pagelist[2]
        fig = pg.figure
        cnv = pg.canvas
        axs = pg.figure.gca()
        fig.clear()
        axs = fig.subplots(2, 1, sharex=True)
        fig.subplots_adjust(right=0.85)
        ax = 1.02
        ay = 1.1
        drag = True
        self.xrr.plot3(axs[0], ax, ay, drag)
        self.xrr.plot4(axs[1], ax, ay, drag)
        # ~ fig.tight_layout()
        cnv.draw()
        
        
    def data_plot1(self):
        pg = self.plotter1.pagelist[0]
        fig = pg.figure
        cnv = pg.canvas
        axs = pg.figure.gca()
        axs.clear()
        fig.subplots_adjust(right=0.85)
        ancx = 1.02
        ancy = 1.1
        drag = True
        self.xrr.plot1(axs, ancx, ancy, drag)
        cnv.draw()
        
        
    def data_plot2(self):
        pg = self.plotter1.pagelist[1]
        fig = pg.figure
        cnv = pg.canvas
        axs = pg.figure.gca()
        
        axs.clear()
        fig.subplots_adjust(right=0.85)
        ancx = 1.02
        ancy = 1.1
        drag = True
        self.xrr.plot2(axs, ancx, ancy, drag)
        cnv.draw()
        
    def data_plot3(self):
        # ~ where()
        cnv = self.plot3.canvas
        fig = self.plot3.figure
        axs = fig.add_subplot(111)
        # ~ axs = fig.gca()
        axs.clear() 
        
        y = self.xrr.ferr    
        n = len(y)
        x = [i for i in range(n)]   
        axs.plot(x, y)
        axs.set_yscale('log')
        fig.tight_layout()
        cnv.draw()

        


    def layer_to_grid(self):
        self.grid_1.ClearGrid()
        nr = self.grid_1.GetNumberRows()
        self.grid_1.DeleteRows( pos=0, numRows=nr)  #, updateLabels=True )
        self.grid_1.SetDefaultCellAlignment(wx.ALIGN_RIGHT, wx.ALIGN_CENTRE)
        self.grid_1.SetRowLabelSize(20)
        
        row = 0
        for L in self.xrr.LL:
            self.grid_1.AppendRows(numRows=1)   #, updateLabels=True)
            self.grid_1.SetRowLabelValue(row, str(row))
            self.grid_1.SetCellValue( row, 0, L._name)
            self.grid_1.SetCellValue( row, 1, L._comp) 
            self.grid_1.SetCellValue( row, 2, "{:.1f}".format(L._dd)) 
            self.grid_1.SetCellValue( row, 3, "{:.1f}".format(L._sg)) 
            self.grid_1.SetCellValue( row, 4, "{:.4f}".format(L._rh)) 
            self.grid_1.SetCellValue( row, 5, "{:.1f}".format(L._Mm)) 
            self.grid_1.SetCellValue( row, 6, "{:.1f}".format(L._Vm))
            row = row +1
        # ~ self.grid_1.AutoSize()
    
    def layer_from_grid(self):
        r = 0
        for L in self.xrr.LL:
            L._name  = self.grid_1.GetCellValue( r, 0)
            L._comp = self.grid_1.GetCellValue( r, 1) 
            L._dd   = float(self.grid_1.GetCellValue( r, 2))
            L._sg   = float(self.grid_1.GetCellValue( r, 3))
            L._rh   = float(self.grid_1.GetCellValue( r, 4))
            r = r + 1
            # ~ print(L._name, L._comp, L._dd, L._sg, L._rh)


            
    def Open_FileDialog(self, ftype = ".xyz"):
        if self.contentNotSaved:
            if wx.MessageBox("Current content has not been saved! Proceed?", "Please confirm",
                             wx.ICON_QUESTION | wx.YES_NO, self) == wx.NO:
                return False
        # otherwise ask the user what new file to open
        # ~ with wx.FileDialog(self, 
        # ~     "Open XYZ file", 
        # ~     wildcard="XYZ files (*.xyz)|*.xyz",
        # ~     style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:
        msg = "Load *"+ftype+" file"
        descript = "Files (*"+ftype+")"
        wildc = descript+"|*"+ftype
        with wx.FileDialog(self, 
            message = msg, 
            wildcard=wildc,
            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fileDialog:
            
            fopen = fileDialog.ShowModal()
            self.pathname    = fileDialog.GetPath()
            self.directory   = fileDialog.GetDirectory()
            self.filename    = fileDialog.GetFilename()
            fileDialog.Destroy()
            if fopen == wx.ID_CANCEL:
                self.pathname  = ""
                self.directory = ""
                self.filename  = ""
                print('!!! cancel open filetype', ftype )
                # ~ check return state in caller event_function and do: event.Skip()
                return False
            else:
                # Proceed loading the file chosen by the user
                return True
                
    def Save_FileDialog(self, ftype = ".xyz"):
        msg = "Save *"+ftype+" file"
        descript = "Files (*"+ftype+")"
        wildc = descript+"|*"+ftype
        with wx.FileDialog(self, 
            message = msg, 
            wildcard=wildc,
            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fileDialog:
                
            fsave = fileDialog.ShowModal()
            self.pathname    = fileDialog.GetPath()
            self.directory   = fileDialog.GetDirectory()
            self.filename    = fileDialog.GetFilename()
            fileDialog.Destroy()
            
            if fsave == wx.ID_CANCEL:
                self.pathname  = ""
                self.directory = ""
                self.filename  = ""
                print('!!! cancel save filetype', ftype )
                return False    # the user changed their mind
            else:
                if ftype not in self.filename:
                    self.filename += ftype
                    self.pathname += ftype
                print(self.filename)
                print(self.pathname)
                return True
    


                
    def Open_try(self):
        try:
            print( "try: open file '%s' " % self.pathname, end='')
            # ~ with open(self.pathname, 'r') as fileh:
                # ~ print(fileh)
            # ~ fileh.close()
            fileh = open(self.pathname,"r")
            fileh.close()
            print( " ... try OK")
            # ~ event.Skip()
            return True
        except IOError:
            print("except: IOError: Cannot open file '%s'." % self.pathname)
            wx.MessageBox("except IOError: Cannot open file '%s' " % self.pathname, 
                "Please confirm",
                wx.ICON_QUESTION | wx.YES_NO, self)
            # ~ event.Skip()
            return False

    def Save_try(self):
            try:
                print( "try: save file :", self.directory, self.filename, end='' )
                file = open(self.pathname,"w")
                file.close()
                print( " ... try OK")
                return True
            except IOError:
                print("Cannot save current data in file '%s'" % self.pathname)
                wx.MessageBox("except IOError: Cannot save file '%s' " % self.pathname, 
                "Please confirm",
                wx.ICON_QUESTION | wx.YES_NO, self)
                return False
    

        
        
    def ev_data_load(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_data_load' not implemented!")
        if not self.Open_FileDialog(ftype = ".dat") :
            event.Skip()
            return
        #else:
        # Proceed loading the file chosen by the user in self.directory
        if self.Open_try():
            self.xrr.data_load(self.pathname)
            self.xrr.data_parse()
            self.xrr.data_init()
            self.data_plot1()
            self.data_plot2()
        event.Skip()

    def ev_layers_load(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_layers_load' not implemented!")
        if not self.Open_FileDialog(ftype = ".xrr") :
            event.Skip()
            return
        #else:
        # Proceed loading the file chosen by the user in self.directory
        if self.Open_try():
            self.xrr.layer_load(self.pathname)
            self.xrr.layer_parse()
            
            self.layer_to_grid()
            self.xrr.layer_init()
            self.xrr.layer_profile()
            self.layer_plot()
            
            # ~ fig = self.plotter2.pagelist[0].figure
            # ~ fig.canvas.draw()
        event.Skip()

    def ev_layers_save(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_layers_save' not implemented!")
        if not self.Save_FileDialog(ftype = ".xrr") :
            print('save file canceled')
            event.Skip()
            return
        #else:
        # Proceed saving the file chosen by the user in self.directory
        if self.Save_try():
            self.xrr.layer_save(self.pathname)
        else:
            print('save try not succesfull')
        event.Skip()

    def ev_layer_add(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_layer_add' not implemented!")
        rows = self.grid_1.GetSelectedRows()
        print('selected rows', rows)
        if len(rows) != 1:
            wx.MessageBox("Select only one row!!!", "Please confirm",
                             wx.ICON_QUESTION | wx.YES_NO, self)
        else:
            r = rows[0]
            self.grid_1.InsertRows(pos=r, numRows=1, updateLabels=True)
            # ~ self.grid_1.SetCellValue( r, 0, 'NoName')
            # ~ self.grid_1.SetCellValue( r, 1, 'Si:1') 
            # ~ self.grid_1.SetCellValue( r, 2, '10.0')   #d 
            # ~ self.grid_1.SetCellValue( r, 3, '1.0')    #sgm 
            # ~ self.grid_1.SetCellValue( r, 4, '2.3')  #rho 
            # ~ self.grid_1.SetCellValue( r, 5, '?')    #Mm
            # ~ self.grid_1.SetCellValue( r, 6, '?')    #Vm

            self.xrr.LL.insert(r, _XLayer('NoName', 'Si;1', d=10.0, s=1.0, r=2.3))
            self.layer_to_grid()

            
        # ~ a = [0,1,2,3,4]
        # ~ a.insert(1, 1.5)
        # ~ b = a[:2] + [3, 3.1, 3.2] + a[2:]
        # ~ c = a[:]
        # ~ c[4:4] = [4.1, 4.2, 4.3]
        # ~ print(a)
        # ~ print(b)
        # ~ print(c)
        event.Skip()

    def ev_layer_del(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_layer_del' not implemented!")
        # ~ self.grid_1.DeleteRows(self, pos=0, numRows=1, updateLabels=True)
        rows = self.grid_1.GetSelectedRows()
        print('selected rows', rows)
        if len(rows) != 1:
            wx.MessageBox("Select only one row!!!", "Please confirm",
                             wx.ICON_QUESTION | wx.YES_NO, self)
        else:
            r = rows[0]
            del self.xrr.LL[r]
            self.layer_to_grid()
            # ~ index = initial_list.index(item1)
            # ~ del initial_list[index]
            # ~ del other_list[index]
        event.Skip()

    def ev_profile_update(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_profile_update' not implemented!")
        self.layer_from_grid()
        self.xrr.layer_init()
        self.layer_to_grid()
        self.xrr.layer_profile()
        self.layer_plot()
        print('xrr()')
        t0 = time.time()
        self.xrr.xrr()
        t1 = time.time()
        print(t1-t0)
        print('plot')
        self.data_plot1()
        self.data_plot2()
        event.Skip()

    def ev_xrr_fit(self, event):  # wxGlade: XRR_Frame.<event_handler>
        # ~ print("Event handler 'ev_xrr_fit' not implemented!")
        self.text_ctrl_log.Clear()

        self.layer_from_grid()
        self.xrr.layer_init()
        self.layer_to_grid()
        
        fitkeys = []
        if self.checkbox_ab.GetValue():            fitkeys.append('ab')
        if self.checkbox_dd.GetValue():            fitkeys.append('dd')
        if self.checkbox_rh.GetValue():            fitkeys.append('rh')
        if self.checkbox_sg.GetValue():            fitkeys.append('sg')
        
        self.xrr.fit_setkeys(fitkeys)
        
        m = self.radio_box_fit.GetSelection()
        if m == 0:  self.xrr.fitmode = 0
        if m == 1:  self.xrr.fitmode = 1
        
        self.xrr.maxfev =100
        self.xrr.fit_xrr()

        self.layer_to_grid()
        self.layer_plot()
        self.data_plot1()
        self.data_plot2()
        event.Skip()
        
        

class MainApp(wx.App):
    """Class Main App."""
    def OnInit(self):
        """Init Main App."""
        self.frame = MainFrame(None, wx.ID_ANY, "")
        self.SetTopWindow(self.frame)
        self.frame.Show(True)
        return True
        
        
        
    

if __name__ == '__main__':
    app = MainApp(0)
    app.MainLoop()

        # ~ self.tx.SetForegroundColour(wx.BLACK)
        
        # ~ font = wx.Font(pointSize=10, 
            # ~ family=wx.MODERN, # DEFAULT, DECORATIVE, ROMAN, SCRIPT, SWISS, MODERN, TELETYPE, MAX
            # ~ style=wx.NORMAL, # FONTSTYLE_ITALIC, FONTSTYLE_SLANT, FONTSTYLE_MAX
            # ~ weight=wx.NORMAL, # THIN, EXTRALIGHT, LIGHT, NORMAL, MEDIUM, SEMIBOLD, BOLD, EXTRABOLD, HEAVY, EXTRAHEAVY, MAX
            # ~ underline=False, 
            # ~ faceName=u'Console')
            
        
         # ~ Font(pointSize, family, style, weight, underline=False,
            # ~ faceName="", encoding=FONTENCODING_DEFAULT)

        # ~ Font(pixelSize, family, style, weight, underline=False,
             # ~ faceName="", encoding=FONTENCODING_DEFAULT)
        
        # Create a font using wx.FontInfo
        

        # ~ FontInfo:
        # ~ AntiAliased     Set anti-aliasing flag.
        # ~ Bold            Use a bold version of the font.
        # ~ Encoding        Set the font encoding to use.
        # ~ FaceName        Set the font face name to use.
        # ~ Family          Set the font family.
        # ~ GetWeightClosestToNumericValue        Get the symbolic weight closest to the given raw weight value.
        # ~ Italic        Use an italic version of the font.
        # ~ Light        Use a lighter version of the font.
        # ~ Slant        Use a slanted version of the font.
        # ~ Strikethrough        Use a strike-through version of the font.
        # ~ Style        Specify the style of the font using one of FontStyle constants.
        # ~ Underlined        Use an underlined version of the font.
        # ~ Weight          Specify the weight of the font.
        
       

        # ~ wx.FONTFAMILY_DEFAULT        Chooses a default font.
        # ~ wx.FONTFAMILY_DECORATIVE    A decorative font.
        # ~ wx.FONTFAMILY_ROMAN         A formal, serif font.
        # ~ wx.FONTFAMILY_SCRIPT        A handwriting font.
        # ~ wx.FONTFAMILY_SWISS         A sans-serif font.
        # ~ wx.FONTFAMILY_MODERN        A fixed pitch font.
        # ~ wx.FONTFAMILY_TELETYPE      A teletype (i.e. monospaced) font.
        # ~ wx.FONTFAMILY_MAX

        # ~ font.setWeighht(wx.FONTWEIGHT_BOLD if cval else wx.FONTWEIGHT_NORMAL)
        # ~ wx.FONTWEIGHT_THIN        Thin font (weight = 100).
        # ~ wx.FONTWEIGHT_EXTRALIGHT    Extra Light (Ultra Light) font (weight = 200).
        # ~ wx.FONTWEIGHT_LIGHT         Light font (weight = 300).
        # ~ wx.FONTWEIGHT_NORMAL    Normal font (weight = 400).
        # ~ wx.FONTWEIGHT_MEDIUM    Medium font (weight = 500).
        # ~ wx.FONTWEIGHT_SEMIBOLD    Semi Bold (Demi Bold) font (weight = 600).
        # ~ wx.FONTWEIGHT_BOLD    Bold font (weight = 700).
        # ~ wx.FONTWEIGHT_EXTRABOLD    Extra Bold (Ultra Bold) font (weight = 800).
        # ~ wx.FONTWEIGHT_HEAVY    Heavy (Black) font (weight = 900).
        # ~ wx.FONTWEIGHT_EXTRAHEAVY    Extra Heavy font (weight = 1000).
        # ~ wx.FONTWEIGHT_MAX


#---------------------------------------------------------------------------
# ~ class PyEvDat(wx.PyEvent):
        # ~ """Simple event to carry arbitrary result data."""
        # ~ def __init__(self, EVT_ID, data):
            # ~ """Init Result Event."""
            # ~ wx.PyEvent.__init__(self)
            # ~ self.SetEventType(EVT_ID)
            # ~ self.data = data
            
            
#===============================================
# ~ https://stackoverflow.com/questions/26312061/what-is-the-difference-between-wx-lib-newevent-newevent-and-wx-neweventtype
# ~ A wx.lib.newevent.NewEvent() is just an easier wxpython way thats been added to make a wx.NewEventType().
# ~ if you have a look at the code in the module newevent you will see what it does.
# ~ """Easy generation of new events classes and binder objects"""
# ~ __author__ = "Miki Tebeka <miki.tebeka@gmail.com>"



# ~ def NewEvent():
    # ~ """Generate new (Event, Binder) tuple
        # ~ e.g. MooEvent, EVT_MOO = NewEvent()
    # ~ """
    # ~ evttype = wx.NewEventType()

    # ~ class _Event(wx.PyEvent):
        # ~ def __init__(self, **kw):
            # ~ wx.PyEvent.__init__(self)
            # ~ self.SetEventType(evttype)
            # ~ self.__dict__.update(kw)

    # ~ return _Event, wx.PyEventBinder(evttype)
    
# ~ https://discuss.wxpython.org/t/wx-pyevent-losing-my-mind/34462/2
# ~ >>> import wx.lib.newevent as ne
# ~ >>> MooEvent, EVT_MOO = ne.NewEvent()
# ~ >>> 
# ~ >>> evt = MooEvent(data=dict(a=1, b=2, c=3))
# ~ >>> 
# ~ >>> evt.data
# ~ {'a': 1, 'b': 2, 'c': 3}
# ~ >>> 
# ~ key = 0
# ~ d = {0:RUN, -1:ABORT, None:END}
# ~ if key in d:
    # ~ func = d[key]
# ~ else:
    # ~ func = RUN
# ~ func()
# ~ d.get(key, RUN)()   # default = RUN
