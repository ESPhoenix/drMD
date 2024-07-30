from typing import Optional


###########################################################################################

def print_drMD_logo() -> None:
    """
    Prints the DRMD logo.

    Returns:
        None
    """

    tealColor = "\033[38;5;37m" 
    resetTextColor = "\033[0m"
    print(tealColor+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
            dddddddd                                                                       
            d::::::d                   MMMMMMMM               MMMMMMMMDDDDDDDDDDDDD         
            d::::::d                   M:::::::M             M:::::::MD::::::::::::DDD     
            d::::::d                   M::::::::M           M::::::::MD:::::::::::::::DD   
            d:::::d                    M:::::::::M         M:::::::::MDDD:::::DDDDD:::::D  
    ddddddddd:::::drrrrr   rrrrrrrrr   M::::::::::M       M::::::::::M  D:::::D    D:::::D 
  dd::::::::::::::dr::::rrr:::::::::r  M:::::::::::M     M:::::::::::M  D:::::D     D:::::D
 d::::::::::::::::dr:::::::::::::::::r M:::::::M::::M   M::::M:::::::M  D:::::D     D:::::D
d:::::::ddddd:::::drr::::::rrrrr::::::rM::::::M M::::M M::::M M::::::M  D:::::D     D:::::D
d::::::d    d:::::d r:::::r     r:::::rM::::::M  M::::M::::M  M::::::M  D:::::D     D:::::D
d:::::d     d:::::d r:::::r     rrrrrrrM::::::M   M:::::::M   M::::::M  D:::::D     D:::::D
d:::::d     d:::::d r:::::r            M::::::M    M:::::M    M::::::M  D:::::D     D:::::D
d:::::d     d:::::d r:::::r            M::::::M     MMMMM     M::::::M  D:::::D    D:::::D 
d::::::ddddd::::::ddr:::::r            M::::::M               M::::::MDDD:::::DDDDD:::::D  
 d:::::::::::::::::dr:::::r            M::::::M               M::::::MD:::::::::::::::DD        
  d:::::::::ddd::::dr:::::r            M::::::M               M::::::MD::::::::::::DDD     
   ddddddddd   dddddrrrrrrr            MMMMMMMM               MMMMMMMMDDDDDDDDDDDDD

⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
                    Molecular Dynamics: Just what the Doctor Ordered!
    """
    +resetTextColor)
###########################################################################################


def print_performing_first_aid() -> None:
    
    yellowText = "\033[33m"
    resetTextColor = "\033[0m"

    print(yellowText+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕

       ██████  ███████ ██████  ███████  ██████  ██████  ███    ███ ██ ███    ██  ██████  
       ██   ██ ██      ██   ██ ██      ██    ██ ██   ██ ████  ████ ██ ████   ██ ██       
       ██████  █████   ██████  █████   ██    ██ ██████  ██ ████ ██ ██ ██ ██  ██ ██   ███ 
       ██      ██      ██   ██ ██      ██    ██ ██   ██ ██  ██  ██ ██ ██  ██ ██ ██    ██ 
       ██      ███████ ██   ██ ██       ██████  ██   ██ ██      ██ ██ ██   ████  ██████  
  
                 ███████ ██ ██████  ███████ ████████      █████  ██ ██████  
                 ██      ██ ██   ██ ██         ██        ██   ██ ██ ██   ██ 
                 █████   ██ ██████  ███████    ██        ███████ ██ ██   ██ 
                 ██      ██ ██   ██      ██    ██        ██   ██ ██ ██   ██ 
                 ██      ██ ██   ██ ███████    ██        ██   ██ ██ ██████  
                                                           
                             Your simulation has crashed...
                             drMD will attempt to rescue it...
         You may need to adjust your simulation parameters in your config file!
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          """
          +resetTextColor)
###########################################################################################
def print_config_error(error: Optional[str] = None) -> None:
    """
    Prints an error message indicating that the config file was not found.

    Returns:
        None
    """
    redText = "\033[31m"
    resetTextColor = "\033[0m"


    print(redText+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          
 ▄████▄  ▒█████   ███▄    █   █████▒██▓ ▄████    ▓█████  ██▀███   ██▀███   ▒█████   ██▀███  
▒██▀ ▀█ ▒██▒  ██▒ ██ ▀█   █ ▓██   ▒▓██▒██▒ ▀█▒   ▓█   ▀ ▓██ ▒ ██▒▓██ ▒ ██▒▒██▒  ██▒▓██ ▒ ██▒
▒▓█    ▄▒██░  ██▒▓██  ▀█ ██▒▒████ ░▒██▒██░▄▄▄░   ▒███   ▓██ ░▄█ ▒▓██ ░▄█ ▒▒██░  ██▒▓██ ░▄█ ▒
▒▓▓▄ ▄██▒██   ██░▓██▒  ▐▌██▒░▓█▒  ░░██░▓█  ██▓   ▒▓█  ▄ ▒██▀▀█▄  ▒██▀▀█▄  ▒██   ██░▒██▀▀█▄  
▒ ▓███▀ ░ ████▓▒░▒██░   ▓██░░▒█░   ░██░▒▓███▀▒   ░▒████▒░██▓ ▒██▒░██▓ ▒██▒░ ████▓▒░░██▓ ▒██▒
░ ░▒ ▒  ░ ▒░▒░▒░ ░ ▒░   ▒ ▒  ▒ ░   ░▓  ░▒   ▒    ░░ ▒░ ░░ ▒▓ ░▒▓░░ ▒▓ ░▒▓░░ ▒░▒░▒░ ░ ▒▓ ░▒▓░
  ░  ▒    ░ ▒ ▒░ ░ ░░   ░ ▒░ ░      ▒ ░ ░   ░     ░ ░  ░  ░▒ ░ ▒░  ░▒ ░ ▒░  ░ ▒ ▒░   ░▒ ░ ▒░
░       ░ ░ ░ ▒     ░   ░ ░  ░ ░    ▒ ░ ░   ░       ░     ░░   ░   ░░   ░ ░ ░ ░ ▒    ░░   ░ 
░ ░         ░ ░           ░         ░       ░       ░  ░   ░        ░         ░ ░     ░     
░          
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          """
          +resetTextColor)
    if error is not None:
      print(f"-->\t{error}")
    exit(1)
###########################################################################################
def print_pdb_error() -> None:
    """
    Prints an error message indicating that problems have been found in the pdb files

    Returns:
        None
    """
    redText = "\033[31m"
    resetTextColor = "\033[0m"


    print(redText+
          """
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕

            ██▓███  ▓█████▄  ▄▄▄▄      ▓█████  ██▀███   ██▀███   ▒█████   ██▀███  
            ▓██░  ██▒▒██▀ ██▌▓█████▄    ▓█   ▀ ▓██ ▒ ██▒▓██ ▒ ██▒▒██▒  ██▒▓██ ▒ ██▒
            ▓██░ ██▓▒░██   █▌▒██▒ ▄██   ▒███   ▓██ ░▄█ ▒▓██ ░▄█ ▒▒██░  ██▒▓██ ░▄█ ▒
            ▒██▄█▓▒ ▒░▓█▄   ▌▒██░█▀     ▒▓█  ▄ ▒██▀▀█▄  ▒██▀▀█▄  ▒██   ██░▒██▀▀█▄  
            ▒██▒ ░  ░░▒████▓ ░▓█  ▀█▓   ░▒████▒░██▓ ▒██▒░██▓ ▒██▒░ ████▓▒░░██▓ ▒██▒
            ▒▓▒░ ░  ░ ▒▒▓  ▒ ░▒▓███▀▒   ░░ ▒░ ░░ ▒▓ ░▒▓░░ ▒▓ ░▒▓░░ ▒░▒░▒░ ░ ▒▓ ░▒▓░
            ░▒ ░      ░ ▒  ▒ ▒░▒   ░     ░ ░  ░  ░▒ ░ ▒░  ░▒ ░ ▒░  ░ ▒ ▒░   ░▒ ░ ▒░
            ░░        ░ ░  ░  ░    ░       ░     ░░   ░   ░░   ░ ░ ░ ░ ▒    ░░   ░ 
                        ░     ░            ░  ░   ░        ░         ░ ░     ░     
⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕⚕
          """
          +resetTextColor)
###########################################################################################



if __name__ == "__main__":
    print_drMD_logo()
    print_config_error()
    print_pdb_error()
    print_performing_first_aid()
  