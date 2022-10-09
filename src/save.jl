#####################################################################
# Data storage
#####################################################################
"stores the presets and supercurrent matrix in two separate CSV files"
function savecsv(p, mat) 
    time_str = string(now())
    pathstr = join([pwd(), "/jcsweepsvsbands/", time_str])
    mkdir(pathstr) 
    CSV.write(join([pathstr,"/presets.csv"]), DataFrame(type2dict(p)))
    CSV.write(join([pathstr,"/jcmat.csv"]), DataFrame(mat, :auto); delim = '\t')
end 


