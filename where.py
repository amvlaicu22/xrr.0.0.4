import sys
def where( msg='', level=1):
    # ~ d = dict( message = msg, 
        # ~ fname  = sys._getframe(1).f_code.co_filename,
        # ~ lineno = sys._getframe(1).f_lineno,
        # ~ func   = sys._getframe(1).f_code.co_name,
        # ~ caller = sys._getframe(2).f_code.co_name)
    # ~ formatstr = "[!]: [{func}] << [{caller}] in [{fname}] : [{lineno}]\n\t >> {message}"
    # ~ textout = formatstr.format( **d)
    # ~ print( textout )
    # ~ or
    fname  = sys._getframe(1).f_code.co_filename
    lineno = sys._getframe(1).f_lineno
    func0  = sys._getframe(1).f_code.co_name
    call1  = sys._getframe(2).f_code.co_name 
    if level == 0:        pass
    if level >= 1:        print(f"{msg}", end=' ')
    if level >= 2:        print(f"in [{func0}]", end=' ')
    if level >= 3:        print(f"called from ({call1})", end=' ')
    if level >= 4:        print(f"in ({fname}) line: ({lineno})", end=' ')
    print()

