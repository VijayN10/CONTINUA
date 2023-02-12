import subprocess
import time

# p = subprocess.Popen([r"E:\Softwares\01\Giraffe\Giraffe.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
# input_file_name = "tex_0"                        # input("Enter the name of the input file: ")
# p.communicate(input=input_file_name.encode())

    
# filename = r"E:\Softwares\01\Giraffe\tex_0\simulation_report.txt"

# while True:
#     with open(filename, "r") as f:
#         lines = f.readlines()
#         if lines:
#             last_line = lines[-1]
#             if "Total solution time:" in last_line:
#                 print("Simulation completed.\n", last_line)
#                 break
#             else:
#                 print(last_line, end="")
#         time.sleep(1)
        
        
        
########################


path = 'E:/Softwares/01/Giraffe/'
folder  = 'E:/Softwares/01/Giraffe/tex_0'
inp  = 'tex_0'

# path of Giraffe.exe
giraffe = path + "Giraffe.exe"

# p = subprocess.Popen(["./Giraffe.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)  # To run from main script 

p = subprocess.Popen([giraffe], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)


input_file_name = inp                                                            # input("Enter the name of the input file: ")
p.communicate(input=input_file_name.encode())

    
# report = r"E:\Softwares\01\Giraffe\tex_0\simulation_report.txt"
report = folder + '/simulation_report.txt'

while True:
    with open(report, "r") as f:
        lines = f.readlines()
        if lines:
            last_line = lines[-1]
            if "Total solution time:" in last_line:
                print("Simulation completed.\n", last_line)
                break
            else:
                print(last_line, end="")
        time.sleep(1) 