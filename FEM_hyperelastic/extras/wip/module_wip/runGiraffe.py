import subprocess
import os

# Run giraffe
#@jit(nopython=False, parallel=True)
def runGiraffe(inp):

    # path : E:/Softwares/01/Giraffe/
    # folder : E:/Softwares/01/Giraffe/tex_0
    # inp  : tex_0

    # path of Giraffe.exe
    # giraffe = path + "Giraffe.exe"
    
    # cur_folder = os.path.dirname(os.getcwd())
    # print(cur_folder)
    
    cur_folder = os.getcwd()
    # print(cur_folder)    
   
    # Changing directory
    os.chdir('Giraffe')

    p = subprocess.Popen(["Giraffe.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    # p = subprocess.Popen([r"E:\Softwares\01\Giraffe\Giraffe.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    input_file_name = inp                                                            # input("Enter the name of the input file: ")
    # p.communicate(input=input_file_name.encode())
    
    output, error = p.communicate(input=input_file_name.encode())
    if error:
        print("Error: ", error.decode())
    else:
        print("Output: ", output.decode())
 
    
    # Going back to current directory or folder
    os.chdir(cur_folder)
    
    opFilePath = 'Giraffe/' + inp + '/monitors/monitor_nodeset_1.txt'

    return opFilePath