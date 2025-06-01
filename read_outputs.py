from pyFAST.input_output.fast_output_file import FASTOutputFile
import os

def read_outfile(c_name, FAST_runDirectory):
    out_file = os.path.join(FAST_runDirectory, c_name + '.out')
    out_sfunc_file = os.path.join(FAST_runDirectory, c_name + 'SFunc.out')
    if os.path.exists(out_file):
        return FASTOutputFile(out_file).toDataFrame()
    elif os.path.exists(out_sfunc_file):
        return FASTOutputFile(out_sfunc_file).toDataFrame()