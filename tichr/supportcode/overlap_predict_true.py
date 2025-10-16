import pandas as pd
import sys

def process_file(filename, n, m): #n和m都是0-based
    df = pd.read_csv(filename,sep="\t",header=None)
    #print(df)
    result = df.groupby(list(df.columns[:n]))[df.columns[m-1]].mean().reset_index()
    return result

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python process_data.py <filename> <n> <m> <outname>")
        sys.exit(1)
    
    #print(sys.argv[2])

    filename = sys.argv[1]
    n = int(sys.argv[2])
    m = int(sys.argv[3])
    outname = str(sys.argv[4])
    

    result = process_file(filename, n, m)
    
    result.to_csv(outname,sep='\t',header=None,index=None)
