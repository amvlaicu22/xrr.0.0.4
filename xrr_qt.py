import sys
import time
import random

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtCore import QThread
from queue import Queue
# ~ from PyQt5.QtCore import Qt, QObject, pyqtSlot
# ~ from PyQt5.QtGui import QIcon, QPalette
# ~ from PyQt5.QtWidgets import (QApplication, QMainWindow,  
# ~ QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout, QFormLayout, 
# ~ QGroupBox, QSizePolicy, QMessageBox, QDialogButtonBox, QLineEdit, 
# ~ QMenu, QMenuBar, QSpinBox, QTextEdit, QFrame, QWidget, QPushButton, 
# ~ QLabel, QComboBox, QDialog, QCheckBox, QRadioButton, QPushButton, 
# ~ QTableWidget, QSlider, QProgressBar)


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

from xrr import XRR, _XLayer
from where import where
from functools import wraps


#########################################################
class StreamStr(QObject):
    signal = pyqtSignal(str)
    def write(self, text): self.newText.emit(str(text))
    def flush(self): pass

class StreamDict(QObject):
    signal = pyqtSignal(dict)
    def write(self, dict_):self.signal.emit(dict_)
    def flush(self): pass

class StreamObj(QObject):
    signal = pyqtSignal(object)
    def write(self, obj):self.signal.emit(obj)
    def flush(self): pass

'''
use in mainframe __init__
    sys.stdout = Stream(signal=self.OnPrint)

    def OnPrint(self, text):
        # self.tx is a QTextEdit
        cursor = self.tx.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text)
        self.tx.setTextCursor(cursor)
        self.tx.ensureCursorVisible()
        # or
        self.tx.insertPlainText( text )
        self.tx.ensureCursorVisible()
'''

# ~ class Qevent_socket(QObject):
    # ~ qsig = pyqtSignal(object)
    # ~ qsig = pyqtSignal(str)
    
    # ~ def __init__(self, qtgui=None, key=1, name='', func=None):
        # ~ QObject.__init__(self)  
        # ~ self.qtgui = qtgui
        # ~ self.key = key
        # ~ self.name = name
        # ~ self.func = func
        # ~ self.data = None
        # ~ self.conn = self.qsig.connect(func)

    # ~ def write(self, data): 
        # ~ self.data = data
        # ~ self.qsig.emit(self)
        # ~ self.qsig.emit(str(data))
    # ~ def flush(self):       pass
    
# ~ class Qevent_socket_dict():
    # ~ def __init__(self, qtgui=None):
        # ~ self.qtgui = qtgui
        # ~ self.dict = {} 
        
    # ~ def add(self, key=0, name='', func=None):
        # ~ self.dict.update( {key : Qevent_socket(self.qtgui, key, name, func) } )
        
    # ~ def post(self, key=0, data={}):
        # ~ self.dict[key].write(data)

class XRR_qt(XRR, QObject):
    signal_run = pyqtSignal(object)
    signal_end = pyqtSignal(object)
    signal_abt = pyqtSignal(object)
    
    # ~ def __init__(self, *args, **kwds):
    def __init__(self):
        print('XRR_wx: __init__()')
        super().__init__()
        # ~ XRR.__init__(self)        
        # ~ QObject.__init__(self)        
        
    def post_event(self, **kw):
        # ~ print(key, data, message, end='', flush=True)
        # ~ # can be redifined when the class is overloaded
        # ~ m = {0:'L-M', 1:'D-E'}
        # ~ d = {'@':'XRR', 'fit':m[self.fitmode], 'fitn':self.fitn}
        # ~ d.update(kw)
        # ~ print(d, end='', flush=True)
        # ~ where(f">>{d}", 2)
        # ~ f ' <text> { <expression> <optional !s, !r, or !a> <optional : format specifier> } <text> ... '
        # ~ Fahrenheit = [32, 60, 102]
        # ~ [f'{((x - 32) * (5/9)):.2f} Celsius' for x in Fahrenheit]
        # ['0.00 Celsius', '15.56 Celsius', '38.89 Celsius']
        d = {'fitn':self.fitn, 'err':f"{self.fiterr:.2e}" }
        d.update(kw)    # add more key:values from **kw
        where(f"\n>>{d}",1)
        
        key = kw['key']
        if key ==  1 : self.signal_run.emit(d)
        if key ==  0 : self.signal_end.emit(d)
        if key == -1 : self.signal_abt.emit(d)
        return
        
# Thread class that executes processing
class QWorkerThread(QThread):
    def __init__(self, gui=None):
        super().__init__()  # ~ QThread.__init__(self) 
        self.qtgui = gui
        print('>>#worker>: __init__()')
        # ~ self.exiting = False
        self.start() # better from the wxgui
        
    def __del__(self):
        # ~ self.exiting = True
        # ~ self.wait()
        print('>>#worker>: __del__()')
                
    # ~ def want_abort(self):
        # ~ self.qtgui.xrr.abort = 1  # just pass it to the xrr
    
    def run(self):
        print('>>#worker>: run()')
        # ~ self.qtgui.xrr.layer_print()
        self.qtgui.xrr.fit()
        
        


# GUI Frame class that spins off the worker thread
class QMainFrame(QMainWindow):
    """Class MainFrame."""
    
    def __init__(self):
        super().__init__()
        self.xrr = XRR_qt()
        self.xrr.signal_run.connect(self.OnEvt_xrr_run)
        self.xrr.signal_end.connect(self.OnEvt_xrr_end)
        self.xrr.signal_abt.connect(self.OnEvt_xrr_abt)

        self.GUI()
        self.Menu()   
        self.statusBar().showMessage('Message in statusbar.')
    
        self.layer_to_grid()
        self.layer_plot()
        self.data_plot1()
        self.data_plot2()
    
        self.show()
        
        self.b1.clicked.connect(self.OnStart)  
        self.b2.clicked.connect(self.OnStop)  
        
        self.worker = None
        
        sys.stdout = StreamObj(signal=self.OnPrintOut)
        sys.stderr = StreamObj(signal=self.OnPrintErr)

    def __del__(self):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
    
    def OnPrintColor(self,r,g,b):
        self.tx.moveCursor(QTextCursor.End)
        self.tx.setTextColor( QColor(r,g,b) )
    
    def OnPrint(self, text):
        self.tx.insertPlainText( text )
        self.tx.ensureCursorVisible()

    def OnPrintOut(self, text):
        font = QFont('Mono', 9, QFont.Light)
        self.tx.setFont(font)
        self.OnPrint(text)

    def OnPrintErr(self, text):
        self.OnPrintColor(255,0,0)
        font = QFont('Mono', 9, QFont.Bold)
        self.tx.setFont(font)
        self.OnPrint(text)
        self.OnPrintColor(0,0,0)
        

    ###################
    
    
    def OnStart(self):
        # Trigger the worker thread unless it's already busy
        if self.worker == None:
            self.tx.clear()
            where("if self.worker == None: Creating worker and start")
            self.xrr.abort = 0
            self.x1.setMinimum(0)
            self.x1.setMaximum(1000)
            self.x1.setValue(0) 
            self.z1.setText('0')
            
            self.layer_from_grid()
            self.xrr.layer_init()
            self.layer_to_grid()
            
            fitkeys = []
            if self.c1.isChecked():            fitkeys.append('ab')
            if self.c2.isChecked():            fitkeys.append('dd')
            if self.c3.isChecked():            fitkeys.append('rh')
            if self.c4.isChecked():            fitkeys.append('sg')

            self.xrr.set_fitkeys(fitkeys)
            
            if self.r1.isChecked():            self.xrr.fitmode = 0
            if self.r2.isChecked():            self.xrr.fitmode = 1
                
            self.b1.setEnabled(False)
            # ~ self.xrr.fit()
            self.worker = QWorkerThread(self)   # self.worker.start()
        else:
            self.OnPrintColor(255,0,0)
            print('!! worker alredy started !!')
            self.OnPrintColor(0,0,0)


            
    def OnStop(self):
        if self.worker == None:
            where("if self.worker == None: enable Fit button.")
            self.b1.setDisabled(False)
            self.b1.setText(' Fit ')
        else:
            if self.worker.isRunning():     # isFinished(), events: started() finished()
                where("if self.worker.isRunning(): Send worker.want_abort()")
                # ~ self.worker.want_abort()
                self.xrr.abort=1
                self.b1.setText(' .... ')
            else:
                where("worker defined but not running: self.worker = None")
                self.b1.setDisabled(False)
                self.b1.setText(' Fit ')
                self.worker = None
               
    
    @pyqtSlot(object)
    def OnEvt_xrr_end(self, sigdict):
        # ~ self.OnPrintColor(155,0,155)
        # ~ where( f"<<{sigdict}", 2 )
        # ~ self.OnPrintColor(0,0,0)
        self.worker = None
        self.b1.setText('Fit')
        self.b1.setDisabled(False)
        self.OnEvt_xrr_update()
        self.xrr.layer_profile()
        self.layer_plot()
    
    
    @pyqtSlot(object)
    def OnEvt_xrr_abt(self, sigdict):
        # ~ self.OnPrintColor(0,0,155)
        # ~ where( f"<<{sigdict}", 2 )
        # ~ self.OnPrintColor(0,0,0)
        pass
    
    
    @pyqtSlot(object)
    def OnEvt_xrr_run(self, sigdict):
        # ~ self.OnPrintColor(0,155,0)
        # ~ where( f"<<{sigdict}", 2)
        # ~ self.OnPrintColor(0,0,0)
        mx = self.x1.maximum()
        nf = sigdict['fitn']
        self.z1.setText(str(nf))
        if (nf  % 100) == 0:
            print('onEvt_xrr_update()')
            self.OnEvt_xrr_update()
        if  nf > mx:
            self.x1.setMaximum(mx+1000)
        self.x1.setValue(nf)

    def OnEvt_xrr_update(self):
        self.layer_to_grid()
        # ~ self.xrr.layer_profile()
        # ~ self.plot_profile()
        self.data_plot1()
        self.data_plot2()
        self.data_plot3()
        # ~ self.p3.axes.clear()
        # ~ self.p3.axes.plot(self.xrr.ferr)
        # ~ self.p3.axes.set_yscale('log')
        # ~ self.p3.draw()
    
    
    def layer_to_grid(self):
        data = []
        data.append(_XLayer._hd)
        for L in self.xrr.LL:
            row = [L._name, L._comp, "{:.1f}".format(L._dd), "{:.1f}".format(L._sg),
            "{:.4f}".format(L._rh), "{:.1f}".format(L._Mm), "{:.1f}".format(L._Vm) ]
            data.append(row)
        self.t1.setData_rows(data)

    
    def layer_from_grid(self):
        for r, L in enumerate(self.xrr.LL):
            L._name  = self.t1.item( r, 0).text()
            L._comp = self.t1.item( r, 1).text() 
            L._dd   = float(self.t1.item( r, 2).text())
            L._sg   = float(self.t1.item( r, 3).text())
            L._rh   = float(self.t1.item( r, 4).text())
            # ~ print(r, L._name, L._comp, L._dd, L._sg, L._rh)
            
    def layer_ins(self):
        row = self.t1.currentRow() 
        self.t1.insertRow(row)
        self.t1.setCurrentCell(row,0)
        self.xrr.LL.insert(row, _XLayer('NoName', 'Si;1', d=10.0, s=1.0, r=2.3))
        self.layer_to_grid()
        
    def layer_del(self):
        row = self.t1.currentRow() 
        self.t1.removeRow(row)
        del self.xrr.LL[row]
        self.layer_to_grid()
    
    def layer_update(self): 
        print('\n>>layer_update(): G2L, L0, L2G, prof, xrr, plot, ', end='') 
        self.layer_from_grid()
        self.xrr.layer_init()
        self.layer_to_grid()
        self.xrr.layer_profile()
        self.xrr.xrr()
        self.layer_plot()
        self.data_plot1()
        self.data_plot2()
        print('<<layer_update')


        
    def layer_plot(self):
        cnv = self.p0
        fig = self.p0.fig
        axs = self.p0.axes        
        fig.clear() 
        axs = fig.subplots(2, 1, sharex=True)
        fig.subplots_adjust(right=0.85)
        ax = 1.02
        ay = 1.1
        drag = True
        self.xrr.plot3(axs[0], ax, ay, drag)
        self.xrr.plot4(axs[1], ax, ay, drag)
        fig.tight_layout()
        cnv.draw()

    def layer_load(self):
        # ~ alert = QMessageBox()
        # ~ alert.setText('layer_load()')
        # ~ alert.exec_()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
            "QFileDialog.getOpenFileName()", "",
            "XRR Files (*.xrr);;All Files (*)", 
            options=options)
        if fileName:
            print(fileName)
            self.xrr.layer_load(fileName)
            self.xrr.layer_parse()
            self.layer_to_grid()
            self.xrr.layer_init()
            # ~ self.xrr.layer_profile()
            # ~ self.layer_plot()
            self.layer_update()
            
        
    def layer_save(self):
        # ~ alert = QMessageBox()
        # ~ alert.setText('layer_save()')
        # ~ alert.exec_()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self,
            "QFileDialog.getSaveFileName()","",
            "XRR Files (*.xrr);;All Files (*)", 
            options=options)
        if fileName:
            print(fileName)
            self.xrr.layer_save(fileName)

    
    def data_load(self):
        # ~ alert = QMessageBox()
        # ~ alert.setText('data_load()')
        # ~ alert.exec_()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self,
            "QFileDialog.getOpenFileName()", "",
            "Data Files (*.dat);;All Files (*)", 
            options=options)
        if fileName:
            print(fileName)
            self.e1.setText(fileName)
            self.xrr.data_load(fileName)
            self.xrr.data_parse()
            self.xrr.data_init()
            self.data_plot1()
            self.data_plot2()

    
    # ~ def openFileNamesDialog(self):
        # ~ options = QFileDialog.Options()
        # ~ options |= QFileDialog.DontUseNativeDialog
        # ~ files, _ = QFileDialog.getOpenFileNames(self,"QFileDialog.getOpenFileNames()", "","All Files (*);;Python Files (*.py)", options=options)
        # ~ if files:
            # ~ print(files)
        
    def data_plot1(self):
        cnv = self.p1
        fig = self.p1.fig
        axs = self.p1.axes
        axs.clear()
        fig.subplots_adjust(right=0.85)
        ancx = 1.02
        ancy = 1.1
        drag = True
        self.xrr.plot1(axs, ancx, ancy, drag)
        cnv.draw()
        
    def data_plot2(self):
        cnv = self.p2
        fig = self.p2.fig
        axs = self.p2.axes
        axs.clear()
        fig.subplots_adjust(right=0.85)
        ancx = 1.02
        ancy = 1.1
        drag = True
        self.xrr.plot2(axs, ancx, ancy, drag)
        cnv.draw()
        
    def data_plot3(self):
        cnv = self.p3
        fig = self.p3.fig
        axs = self.p3.axes
        axs.clear() 
        y = self.xrr.ferr    
        n = len(y)
        x = [i for i in range(n)]   
        axs.plot(x, y)
        axs.set_yscale('log')
        fig.tight_layout()
        cnv.draw()
    

    def QB(self, label, tip='tooltip', hsize=140, vsize=25, func=None): #, signal='clicked()'):
        button = QPushButton(label)
        button.setToolTip(tip)
        button.resize(hsize,vsize)
        if func is not None:
            button.clicked.connect(func)
        return button
        
    def GUI(self):
        self.setGeometry(80, 25, 1280, 730)
        # ~ self.setWindowTitle("my first window")
        self.setWindowTitle("X-ray reflectivity - Parrat algorithm")

        
        L = QGridLayout()
        F = QFrame(self)
        self.setCentralWidget(F)
        F.setLayout(L)

        L0_0 = QHBoxLayout()
        L0_0.addWidget(QLabel('Data') )
        L0_0.addWidget(self.QB('load', func = self.data_load) )
        self.e1 = QLineEdit();           
        L0_0.addWidget(self.e1)
        F0_0 = QFrame()
        F0_0.setLayout(L0_0)


        L0_1 = QHBoxLayout()
        L0_1.addWidget(QLabel('Multi-Layers'))
        L0_1.addWidget(self.QB( 'load', func = self.layer_load))
        L0_1.addWidget(self.QB( 'save', func = self.layer_save))
        F0_1 = QFrame()
        F0_1.setLayout(L0_1)
        
        L1_0 = QVBoxLayout()
        self.tb = MyTabWidget()
        L1_0.addWidget(self.tb)
        F1_0 = QFrame()
        F1_0.setLayout(L1_0)
        self.p1 = self.tb.cnv1
        self.p2 = self.tb.cnv2
        self.p0 = self.tb.cnv3


        L1_1 = QVBoxLayout()
        self.t1 = MyTableWidget(self)
        L1_1.addWidget(self.t1)
        F1_1 = QFrame()
        F1_1.setLayout(L1_1)


        L2_0 = QHBoxLayout()
        self.b1 = self.QB( 'Fit', func = self.OnStart)
        self.b2 = self.QB( 'stop', func = self.OnStop)
        L2_0.addWidget(self.b1)
        L2_0.addWidget(self.b2)
        self.r1 = QRadioButton('L-M');  self.r1.setChecked(True);   
        self.r2 = QRadioButton('D-E');  self.r2.setChecked(False); 
        self.z1 = self.QB('0');    self.z1.resize(100, 20)
        # ~ self.x1 = QProgressBar();       self.x1.resize(120, 20)
        self.x1 = QSlider(Qt.Horizontal);       self.x1.resize(120, 20)
        self.x1.setMinimum(0)
        self.x1.setMaximum(1000)
        self.x1.setValue(0)
        self.x1.setEnabled(False)
        self.x1.setTickPosition(QSlider.TicksBothSides) # TicksBelow TicksLeft
      
        self.c1 = QCheckBox('a0b0');    self.c1.setChecked(True)
        self.c2 = QCheckBox('dd');      self.c2.setChecked(True)
        self.c3 = QCheckBox('sg');      self.c3.setChecked(True)
        self.c4 = QCheckBox('rh');      self.c4.setChecked(True)
        L2_0.addWidget(self.r1)  
        L2_0.addWidget(self.r2)
        L2_0.addWidget(self.z1)
        L2_0.addWidget(self.x1)
        L2_0.addWidget(self.c1)
        L2_0.addWidget(self.c2)
        L2_0.addWidget(self.c3)
        L2_0.addWidget(self.c4)
        F2_0 = QFrame()
        F2_0.setLayout(L2_0)
        



        L2_1 = QHBoxLayout()
        L2_1.addWidget(self.QB( 'insert', func = self.layer_ins) )
        L2_1.addWidget(self.QB( 'delete', func = self.layer_del) )
        L2_1.addWidget(self.QB( 'update', func = self.layer_update) )
        F2_1 = QFrame()
        F2_1.setLayout(L2_1)
        
        L3_1 = QHBoxLayout()
        self.p3 = PlotCanvas(self, title='fitting error', width=3, height=1)
        self.p3.fig.tight_layout()

        L3_1.addWidget(self.p3)
        F3_1 = QFrame()
        F3_1.setLayout(L3_1)

        L4_1 = QHBoxLayout()
        self.tx = QTextEdit()       
        L4_1.addWidget(self.tx)
        F4_1 = QFrame()
        F4_1.setLayout(L4_1)
      
        L.addWidget(F0_0, *(0, 0))
        L.addWidget(F2_0, *(1, 0))             
        L.addWidget(F1_0, *(2, 0))
        L.addWidget(F3_1, *(3, 0))
            
        L.addWidget(F0_1, *(0, 1))
        L.addWidget(F2_1, *(1, 1))
        L.addWidget(F1_1, *(2, 1))
        L.addWidget(F4_1, *(3, 1))
        
        
    def Menu(self):
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('File')
        editMenu = mainMenu.addMenu('Edit')
        viewMenu = mainMenu.addMenu('View')
        searchMenu = mainMenu.addMenu('Search')
        toolsMenu = mainMenu.addMenu('Tools')
        helpMenu = mainMenu.addMenu('Help')
        
        exitButton = QAction(QIcon('exit24.png'), 'Exit', self)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(self.close)
        fileMenu.addAction(exitButton)
        
        # ~ self.show()
        
    

class PlotCanvas(FigureCanvas):

    def __init__(self, parent=None, width=4, height=3, dpi=100, title='PlotCanvas'):

        self.title = title
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.axes.set_title(self.title)

        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,    QSizePolicy.Expanding,      QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
        



# ~ from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
# ~ from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
# ~ from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

class MyTabWidget(QWidget):
    
    # ~ def __init__(self, parent):
        # ~ super(QWidget, self).__init__(parent)
    def __init__(self):
        super(QWidget, self).__init__()
        
        # ~ # Initialize tab screen
        self.tabs = QTabWidget()
        self.tabs.resize(300,200)
        self.tabs.setTabPosition(2)
        
        # ~ # Add tabs to widget
        self.layout = QVBoxLayout(self)
        self.layout.addWidget(self.tabs)
        self.setLayout(self.layout)
        
        # ~ # Add tabs
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tabs.addTab(self.tab1,"I vs 2theta")
        self.tabs.addTab(self.tab2,"I*Q^4 vs Q")
        self.tabs.addTab(self.tab3,"Layer profile")
        
        # ~ # Create first tab
        self.cnv1 = PlotCanvas(self, title='c1')
        self.nav1 = NavigationToolbar(self.cnv1, self)
        self.lay1 = QVBoxLayout(self)
        self.lay1.addWidget(self.cnv1)
        self.lay1.addWidget(self.nav1)
        self.tab1.setLayout(self.lay1)
        
        # ~ # Create second tab
        self.cnv2 = PlotCanvas(self, title='c2')
        self.nav2 = NavigationToolbar(self.cnv2, self) 
        self.lay2 = QVBoxLayout(self)
        self.lay2.addWidget(self.cnv2)
        self.lay2.addWidget(self.nav2)
        self.tab2.setLayout(self.lay2)

        # ~ # Create third tab
        self.cnv3 = PlotCanvas(self, title='c2')
        self.nav3 = NavigationToolbar(self.cnv3, self) 
        self.lay3 = QVBoxLayout(self)
        self.lay3.addWidget(self.cnv3)
        self.lay3.addWidget(self.nav3)
        self.tab3.setLayout(self.lay3)
        

class MyTableWidget(QTableWidget):
    data_row = []
    data_col = []
            
    def __init__(self, data, *args):
        QTableWidget.__init__(self, *args)
        self.resize(300,200)
        
        self.setShowGrid(True)
        self.setRowCount(0)
        self.setColumnCount(0)
        self.setData_rows()

    def setData_rows(self, 
        data=[[   'name',     'comp.',        'd [A]',    's [A]',    'g/cm^3',   'g/mol',    'A^3'],
              [   'glass',    'Si;1;O;2',    '0',         '5.0',      '2.6',      '',         ''],
              [   'W',        'W;1',          '10.0',     '5.0',      '20.0',     '',         '']]): 
        self.data_row = data
        nr = len(self.data_row) -1 # 1st row is for header
        nc = len(self.data_row[0])
        self.setRowCount(nr)
        self.setColumnCount(nc)
        self.setHorizontalHeaderLabels(self.data_row[0])
        for r, row in enumerate(self.data_row[1:], start=0):
            for c, item in enumerate(row):
                newitem = QTableWidgetItem(item)
                self.setItem(r, c, newitem)
        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        
    @pyqtSlot()
    def on_click(self):
        print("\n")
        for item in self.tableWidget.selectedItems():
            print(item.row(), item.column(), item.text())
        
    def set_attrib(self):
        self.horizontalHeader().setDefaultSectionSize(80)
        self.horizontalHeader().resizeSection(0, 100)
        self.verticalHeader().setVisible(True)
        self.verticalHeader().setDefaultSectionSize(19)
        self.horizontalHeader().setSectionResizeMode(0, QHeaderView.Interactive)
        self.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.resizeColumnsToContents()
        self.resizeRowsToContents()
        self.setItem(1,2, QTableWidgetItem("Table Cell"))
    
    def setData_cols(self, 
        data={'col1':['1','2','3','4'],
        'col2':['1','2','1','3'],
        'col3':['1','1','2','1']}): 
        self.data_col = data
        horHeaders = []
        for n, key in enumerate(sorted(self.data_col.keys())):
            horHeaders.append(key)
            for m, item in enumerate(self.data[key]):
                newitem = QTableWidgetItem(item)
                self.setItem(m, n, newitem)
        self.setHorizontalHeaderLabels(horHeaders)



if __name__ == '__main__':

    app = QApplication(sys.argv)
    # ~ app = QApplication([])
    # ~ app.setStyle('Fusion')
    # ~ 'Fusion' , 'Windows' , 'WindowsVista' Windows only and 'Macintosh' Mac only.
    
    # ~ app.setStyleSheet("QPushButton { margin: 1ex; }")
    
    # ~ palette = QPalette()
    # ~ palette.setColor(QPalette.ButtonText, Qt.blue)
    # ~ app.setPalette(palette)
    
    frame = QMainFrame()
    sys.exit(app.exec_())
    


"""
# ===============================================================
# The new Stream Object which replaces the default stream 
# associated with sys.stdout
# This object just puts data in a queue!

class WriteStream(object):
    def __init__(self,queue):   self.queue = queue
    def write(self, text):      self.queue.put(text)
    def flush(self):            pass;
        
# A QObject (to be run in a QThread) which sits waiting for data to come through a Queue.Queue().
# It blocks until data is available, and one it has got something from the queue, it sends
# it to the "MainThread" by emitting a Qt Signal 

class MyReceiver(QObject):
    mysignal = pyqtSignal(str)
    def __init__(self,queue,*args,**kwargs):
        QObject.__init__(self,*args,**kwargs)
        self.queue = queue
        
    @pyqtSlot()
    def run(self):
        while True:
            text = self.queue.get()
            self.mysignal.emit(text)
            
'''
    # Within your main window class...
    def __init__(self, parent=None, **kwargs):
        # ...
        # Install the custom output stream
        self.redirect()
        
    def __del__(self):
        # Restore sys.stdout
        sys.stdout = sys.__stdout__

    def redirect(self):
        # Create Queue and redirect sys.stdout to this queue
        self.queue = Queue()
        # Create thread that will listen on the other end of the queue, 
        # ~ and send the text to the textedit in our application
        self.thread = QThread()
        self.my_receiver = MyReceiver(self.queue)
        self.my_receiver.mysignal.connect(self.append_text)
        self.my_receiver.moveToThread(self.thread)
        self.thread.started.connect(self.my_receiver.run)
        self.thread.start()
        sys.stdout = WriteStream(self.queue)
        sys.stderr = WriteStream(self.queue)
'''
# ===============================================================

def _get_format_from_style(self, token, style):
        ''' Returns a QTextCharFormat for token by reading a Pygments style.
        '''
        result = QtGui.QTextCharFormat()
        for key, value in list(style.style_for_token(token).items()):
            if value:
                if key == 'color':
                    result.setForeground(self._get_brush(value))
                elif key == 'bgcolor':
                    result.setBackground(self._get_brush(value))
                elif key == 'bold':
                    result.setFontWeight(QtGui.QFont.Bold)
                elif key == 'italic':
                    result.setFontItalic(True)
                elif key == 'underline':
                    result.setUnderlineStyle(
                        QtGui.QTextCharFormat.SingleUnderline)
                elif key == 'sans':
                    result.setFontStyleHint(QtGui.QFont.SansSerif)
                elif key == 'roman':
                    result.setFontStyleHint(QtGui.QFont.Times)
                elif key == 'mono':
                    result.setFontStyleHint(QtGui.QFont.TypeWriter)
        return result 
        
# ~ class StdoutWrapper(object):
# ~ class StdoutWrapper(QTextEdit):
    # ~ def __init__(self, outwidget):
        # ~ super().__init__() # if using QTextEdit
        # ~ self.widget = outwidget
    # ~ def write(self, s):
        # ~ self.widget.moveCursor(QTextCursor.End)
        # ~ self.widget.insertText( s )
    # ~ def flush(self):
        # ~ pass;
# ~ use:
# ~ sys.stdout = StdoutWrapper(self.tx)
# ~ sys.stderr = StdoutWrapper(self.tx)
# similar for stderr, but you might want an error dialog or make 
# the text stand out using appendHtml
        
'''
import sys
orig_stdout = sys.stdout
f = file('out.txt', 'w')
sys.stdout = f
for i in range(2):
    print ('i = ', i)
sys.stdout = orig_stdout
f.close()
'''

class timer_class():
    def __init__():
        self.timer = QBasicTimer()
        self.step = 0
        
    def timerEvent(self, e):
        if self.step >= 100:
            # ~ self.timer.stop()
            # ~ self.btn.setText('Finished')
            self.step = 0
            return
            
        self.step = self.step + 1
        self.x1.setValue(self.step)
        
    def timeTest():
        if self.timer.isActive():
            self.timer.stop()
            self.step = 0
            self.b1.setText('Fit')
        else:
            self.timer.start(100, self)
            self.b1.setText('Stop')
        alert = QMessageBox()
        alert.setText('data_fit()')
        alert.exec_()
        return
    
        
"""

'''
# https://www.riverbankcomputing.com/static/Docs/PyQt5/signals_slots.html

class on_write_qt_emit(QObject):
    textWritten = pyqtSignal(str, int)
    def __init__(self, cbfunc, cval, *args,**kwargs):
        QObject.__init__(self,*args,**kwargs)
        self.textWritten.connect(cbfunc)
        self.cval = cval # 0:black, 1: red
    def write(self, text):
        self.textWritten.emit(str(text), self.cval)
    def flush(self):       pass;
'''

'''
# https://www.learnpyqt.com/examples/create-desktop-weather-app/
class WorkerSignals(QObject):

    #Defines the signals available from a running worker thread.
 
    finished = pyqtSignal()
    error = pyqtSignal(str)
    result = pyqtSignal(dict, dict)
    
class WeatherWorker(QRunnable):

    #Worker thread for weather updates.

    signals = WorkerSignals()
    is_interrupted = False

    def __init__(self, location):
        super(WeatherWorker, self).__init__()
        self.location = location

    @pyqtSlot()
    def run(self):
        try:
            params = dict(
                q=self.location,
                appid=OPENWEATHERMAP_API_KEY
            )

            url = 'http://api.openweathermap.org/data/2.5/weather?%s&units=metric' % urlencode(params)
            r = requests.get(url)
            weather = json.loads(r.text)

            # Check if we had a failure (the forecast will fail in the same way).
            if weather['cod'] != 200:
                raise Exception(weather['message'])

            url = 'http://api.openweathermap.org/data/2.5/forecast?%s&units=metric' % urlencode(params)
            r = requests.get(url)
            forecast = json.loads(r.text)

            self.signals.result.emit(weather, forecast)

        except Exception as e:
            self.signals.error.emit(str(e))

        self.signals.finished.emit()
        
class MainWindow(QMainWindow, Ui_MainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.pushButton.pressed.connect(self.update_weather)
        self.threadpool = QThreadPool()
        self.show()
        
    def update_weather(self):
        worker = WeatherWorker(self.lineEdit.text())
        worker.signals.result.connect(self.weather_result)
        worker.signals.error.connect(self.alert)
        self.threadpool.start(worker)

    def alert(self, message):
        alert = QMessageBox.warning(self, "Warning", message)
        
    def weather_result(self, weather, forecasts):
        self.latitudeLabel.setText("%.2f °" % weather['coord']['lat'])
        self.longitudeLabel.setText("%.2f °" % weather['coord']['lon'])

        self.windLabel.setText("%.2f m/s" % weather['wind']['speed'])

        self.temperatureLabel.setText("%.1f °C" % weather['main']['temp'])
        self.pressureLabel.setText("%d" % weather['main']['pressure'])
        self.humidityLabel.setText("%d" % weather['main']['humidity'])

        self.weatherLabel.setText("%s (%s)" % (
            weather['weather'][0]['main'],
            weather['weather'][0]['description']
        )

'''

'''
# https://stackoverflow.com/questions/43964766/pyqt-emit-signal-with-dict
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QApplication


class Emiterer(QtCore.QThread):
    f = QtCore.pyqtSignal(dict)

    def __init__(self):
        super(Emiterer, self).__init__()

    def run(self):
        self.f.emit({"2": {}})
        # self.f.emit({2: {}})  < == this don't work!


class Main(QtWidgets.QMainWindow):

    def __init__(self):
        super(Main, self).__init__()
        self.e = Emiterer()
        self.e.f.connect(self.finised)
        self.e.start()

    def finised(self, r_dict):
        print(r_dict)


if __name__ == '__main__':
    import sys
    app = QApplication(sys.argv)
    m = Main()
    m.show()
    sys.exit(app.exec_())
    


Use object instead of dict in the pyqtSignal definition. E.g.

class Emiterer(QtCore.QThread):
    f = QtCore.pyqtSignal(object)

The reason for this is that signals defined as pyqtSignal(dict) are actually interpreted the same as pyqtSignal('QVariantMap') by PyQt5 and QVariantMap can only have strings as keys.

You can check this (for your specific class) with

m = Emiterer.staticMetaObject
method = m.method(m.methodCount() - 1)  # the last defined signal or slot
print(method.methodSignature())

This would print PyQt5.QtCore.QByteArray(b'f(QVariantMap)')

'''

'''

class on_write_qt_emit(QObject):
    textWritten = pyqtSignal(str, int)
    def __init__(self, cbfunc, cval, *args,**kwargs):
        QObject.__init__(self,*args,**kwargs)
        self.textWritten.connect(cbfunc)
        self.cval = cval # 0:black, 1: red
    def write(self, text):
        # self.textWritten.emit(str(text), self.cval)
    def flush(self):       pass;
    

    # Within your main window class...
    def __init__(self, parent=None, **kwargs):
        # ...
        # Install the custom output stream
        sys.stdout = on_write_qt_emit(textWritten=self.append2QTextEdit)

    def __del__(self):
        # Restore sys.stdout
        sys.stdout = sys.__stdout__

    def append2QTextEdit(self, text):
        """Append text to the QTextEdit."""
        # Maybe QTextEdit.append() works as well, but this is how I do it:
        cursor = self.textEdit.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.textEdit.setTextCursor(cursor)
        self.textEdit.ensureCursorVisible()
'''

'''
print('___________________________________________')
import sys
import inspect
# https://www.oreilly.com/library/view/python-cookbook/0596001673/ch14s08.html
# ~ sys._getframe. 
# ~ This function returns a frame object whose attribute f_code is a code object 
# ~ and the co_name attribute of that object is the function name

this_function_name = sys._getframe(  ).f_code.co_name
this_line_number = sys._getframe(  ).f_lineno
this_filename = sys._getframe(  ).f_code.co_filename

def whoami(  ):
    import sys
    return sys._getframe(1).f_code.co_name
def callersname(  ):
    import sys
    return sys._getframe(2).f_code.co_name

def where():
    import sys
    print( sys._getframe(1).f_code.co_filename, ':',
        sys._getframe(1).f_lineno, ':', 
        sys._getframe(1).f_code.co_name, '<<',
        sys._getframe(2).f_code.co_name)
    
def foo():
    where()
    
foo()


# ~ https://medium.com/@vadimpushtaev/name-of-python-function-e6d650806c4

def fname(func):
    @wraps(func)
    def tmp(*args, **kwargs):
        print('>>', func.__name__ , '<<')
        return func(*args, **kwargs)
    return tmp

fname1 = lambda n=0: sys._getframe(n + 1).f_code.co_name
# for current func name, specify 0 or no argument.
# for name of caller of current func, specify 1.
# for name of caller of caller of current func, specify 2. etc.

fname2 = lambda: inspect.stack()[1][3]

def fname3(): 
    print(inspect.stack()[0][3])
    print(inspect.stack()[1][3]) #will give the caller of foos name, if something called foo
    
    print(inspect.stack()[0].function, 'inspect.stack[0].function')
    print(inspect.stack()[1].function, 'inspect.stack[1].function') 
    
    print( fname1(),  'currentFuncName() with sys._getframe(n)')
    print( fname1(1), 'currentFuncName(1)' )   
    
    print(inspect.currentframe().f_code.co_name, 'inspect.currentframe().f_code.co_name')
    print(inspect.currentframe().f_back.f_code.co_name, 'inspect.currentframe().f_back.f_code.co_name')
    # ~ print( "my name is '{}'".format( variable )

@fname
def my_funky_name():
    print( "STUB" )
    fname1()
    fname2()
    fname3()

my_funky_name()

# ~ $ python -m timeit -s 'import inspect, sys' 'inspect.stack()[0][0].f_code.co_name'
# ~ 1000 loops, best of 3: 499 usec per loop
# ~ $ python -m timeit -s 'import inspect, sys' 'inspect.stack()[0][3]'
# ~ 1000 loops, best of 3: 497 usec per loop
# ~ $ python -m timeit -s 'import inspect, sys' 'inspect.currentframe().f_code.co_name'
# ~ 10000000 loops, best of 3: 0.1 usec per loop
# ~ $ python -m timeit -s 'import inspect, sys' 'sys._getframe().f_code.co_name'
# ~ 10000000 loops, best of 3: 0.135 usec per loop

def safe_super(_class, _inst):
    """safe super call"""
    try:
        return getattr(super(_class, _inst), _inst.__fname__)
    except:
        return (lambda *x,**kx: None)


def with_name(function):
    def wrap(self, *args, **kwargs):
        self.__fname__ = function.__name__
        return function(self, *args, **kwargs)
    return wrap

class A(object):

    def __init__():
        super(A, self).__init__()

    @with_name
    def test(self):
        print 'called from A\n'
        safe_super(A, self)()

class B(object):

    def __init__():
        super(B, self).__init__()

    @with_name
    def test(self):
        print 'called from B\n'
        safe_super(B, self)()

class C(A, B):

    def __init__():
        super(C, self).__init__()

    @with_name
    def test(self):
        print 'called from C\n'
        safe_super(C, self)()

# ~ testing it :

a = C()
a.test()

# ~ Inside each @with_name decorated method you have access to 
# ~ self.__fname__ as the current function name.

def foo():
    """foo docstring"""
    print(eval(sys._getframe().f_code.co_name).__doc__)
    
foo()


# ~ print('___________________________________________')
'''

        # ~ radiobutton = QRadioButton("Australia")
        # ~ radiobutton.setChecked(True)
        # ~ radiobutton.country = "Australia"
        # ~ radiobutton.toggled.connect(self.onClicked)
        # ~ self.show()
