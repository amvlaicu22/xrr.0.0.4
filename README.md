# XRR
X-ray reflectivity modeling for thin film multi-layer structures using Parrat algorithm 

### Features
* Layer model fitting
* Density Profile fitting
* Fitting algorithm:
	* [Levenberg-Marquardt](
	https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm)
	* [Differential-Evolution](
	https://en.wikipedia.org/wiki/Differential_evolution)

Future developments:
* atomic composition profile fitting
* console argument parsing
	
#### License

This project is licensed under the **GNU General Public License version 2** 

* [https://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html](
https://www.gnu.org/licenses/old-licenses/gpl-2.0-standalone.html)
* [https://opensource.org/licenses/GPL-2.0](
https://opensource.org/licenses/GPL-2.0)
* see the local [LICENSE](./LICENSE) file for details

#### Acknowledgments
* Cromer-Liberman code for Igor Pro ver. 4.0 and higher, translated from original Fortran code by **Jan Ilavsky**, 9/18/2006:
	* [https://usaxs.xray.aps.anl.gov/software/cromer-liberman](
	https://usaxs.xray.aps.anl.gov/software/cromer-liberman)


#### Prerequisites and Installing

* [python](https://www.python.org/) 3.6 or higher
* package manager [pip](https://pip.pypa.io/en/stable/)
* [numpy](https://pypi.org/project/numpy/)
* [scipy](https://pypi.org/project/scipy/)
* [mpmath](https://pypi.org/project/mpmath/) for the **erf** function
* [matplotlib](https://pypi.org/project/matplotlib/)
```
sudo -H pip3 install numpy scipy mpmath matplotlib
```
For Graphical User Interface
* [wxPython](https://pypi.org/project/wxPython/), [PyQt5](https://pypi.org/project/PyQt5/) 
```
sudo -H pip3 install wxPython PyQt5
```

##### Install from PyPi
[https://test.pypi.org/project/xrr-amvlaicu22/](
https://test.pypi.org/project/xrr-amvlaicu22/) 

try it in your local folder before using sudo -H:

```bash
pip3 install -i https://test.pypi.org/simple/ --no-deps xrr-amvlaicu22
```
##### Download from github
[https://github.com/amvlaicu22/xrr](
https://github.com/amvlaicu22/xrr)
##### Files
* xrr.py 	: console program
* CromerLiberman.py : atomic scattering factors 
* xrr_wx2.py 	: wx interface 
* xrr_wx2gui.py	: wxglade generated GUI 
* xrr_wx2gui.wxg : wxglade file for GUI 
* xrr_qt2.py 	: Qt5 interface 
* xrr_qt2gui.ui : qt-designer GUI file
 

##### Running
```bash
python3 xrr.py
python3 xrr_wx2.py
python3 xrr_qt2.py
    
```
##### *test data*:
* xrr_test.dat	: xrr data (2tht, int)
* xrr_test.xrr	: layer model

#### Author(s) (contact)
* *Vlaicu Aurel-Mihai* [(amvlaicu22 at yahoo.com)](mailto:amvlaicu22_at_yahoo.com)

