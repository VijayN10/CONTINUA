import os

# Check the size of output file

# The function returns a boolean value indicating whether the file is empty or not. 
# If the file is empty, os.path.getsize('.txt') will return 0 and checkGiraffeOutputs will return True. 
# If the file is non-empty, os.path.getsize('.txt') will return a value greater than 0 and checkGiraffeOutputs will return False.
#@jit(nopython=False, parallel=True)
def checkGiraffeOutputs(opFilePath):

    if os.path.getsize(opFilePath) == 0:
        print("File is empty! We are running the Giraffe again...")
        return True
    
    return False