
import numpy as np
import scipy
from uonidtoolbox import _demos

def demo(demo_type=[]):
    if not demo_type:
        demo_number = -2
        while demo_number == -2:
            demo_number = _requestDemoFromUser() # waits for user to select
            
            if demo_number == -1: # then quit was requested
                return

            if not (0 <= demo_number < len(_demo_map)):
                print("-"*50+"\nNot a valid option. Pick from the list")
                demo_number = -2
        #endwhile

        print(f"\nRunning '{_demo_map[demo_number][1]}' demo...")
        demo(_demo_map[demo_number][0])
    else:
        if isinstance(demo_type, int):
            if 0 < demo_type < len(_demo_map):
                demo_type = _demo_map[demo_number][0]
            else:
                raise Exception("not a valid option, call this function with no argument for help")

        if not isinstance(demo_type, str):
            raise Exception("not a valid option, call this function with no argument for help")

        if demo_type not in [t[0] for t in _demo_map]:
            raise Exception("not a valid option, call this function with no argument for help")

        getattr(_demos, 'demo_'+demo_type)()
#endfunction


def _requestDemoFromUser():
    _printUONAscii()
    print("\nUoN ID Toolbox: List of demos\n" + "="*50)
    for i in range(1, len(_demo_map)):
        print(f"{i:<3}: "+_demo_map[i][1])

    userin = input("\nPlease select a demo from the list by its index (type 'q' to quit)': ")
    
    # handle empty input
    if len(userin) < 1: raise Exception("-"*50+"\nPick a number buddy, or type 'q' to quit\n")

    # handle quit request
    if userin[0].lower() == 'q': return -1

    # see if the input can be converted to type int
    try:
        num = int(userin)
    except: 
        print("-"*50+"\nNeeds to be a number buddy, or type 'q' to quit")
        return -2

    if num < 0:
        print("-"*50+"\nPick a number from the list, or type 'q' to quit")
        return -2

    return num
#endfunction

_demo_map = [
('', "empty"), 
('ar', "Autoregressive"), 
('arx', "Autoregressive with exogenous input"), 
('fir', "Finite impulse response"),
('oe', "Output-error"),
('bj', "Box-Jenkins"),
]

def _printUONAscii():
    print("")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣀⣤⡶⠒⢀⣠⣴⣶⣿⣿⣿⠟⠋⠉⠀⠀⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⠀⠀⠀⠀⢀⣠⣶⣿⡿⠉⣠⣶⣿⣿⣿⣿⣿⣿⣿⣆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⠀⠀⢀⣴⣿⣿⣿⠟⣠⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⠀⣰⣿⣿⣿⣿⣿⣾⣿⣿⣿⣿⣿⣿⡿⠛⠉⠉⣉⣙⣿⣿⣷⡄⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⣼⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠋⠀⠀⣰⣾⣿⣿⣿⣿⣿⣿⣦⡀⠀⠀⠀⠀")
    print("⠀⣸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠇⠀⠀⠀⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣄⠀⠀⠀")
    print("⢀⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡆⠀⠀⠀⠘⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣄⠀")
    print("⢸⣿⣿⣿⣿⣿⣿⣿⡿⠟⢛⣻⣿⣿⣷⡀⠀⠀⠀⠀⠈⠉⠉⠛⠛⠻⢿⣿⣿⣟⠿⣷")
    print("⢸⣿⣿⣿⣿⣿⣿⠏⠀⢠⣿⣿⣿⣿⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⣿⠿⠃⠀")
    print("⠘⣿⣿⣿⣿⣿⡏⠀⠀⢸⣿⣿⣿⣿⣿⣿⣿⣷⣤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀")
    print("⠀⢻⣿⣿⣿⣿⡇⠀⠀⠀⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣦⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀")
    print("⠀⠈⢿⣿⣿⣿⣿⡄⠀⠀⠀⠈⠙⠻⢿⣿⣿⣿⣿⣿⣿⣿⣿⣷⡄⠀⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⠈⢿⣿⣿⣿⣿⣆⠀⠀⠀⠀⠀⠀⠀⠉⠙⠻⢿⣿⣿⣿⣿⣿⡄⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⠀⠀⠹⣿⣿⣿⣿⣷⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⣿⣿⣿⡇⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⠀⠀⠀⠈⠛⢿⣿⣿⣿⣷⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⣿⡿⠁⠀⠀⠀⠀⠀⠀")
    print("⠀⠀⠀⠀⠀⠀⠀⠀⠈⠛⠻⠿⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠼⠋⠀⠀⠀⠀⠀⠀⠀⠀")
