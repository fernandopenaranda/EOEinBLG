"stores the presets, sweep vector and spectrum matrix as three CSV files stored 
in a new folder whose name is the output of `string(now())`"
function savecsv(p, mat, excitation_energy) 
    time_str = string(now())
    pathstr = join([pwd(), "/jcsweepsvsbands/", time_str])
    mkdir(pathstr) 
    CSV.write(join([pathstr,"/presets.csv"]), DataFrame(type2dict(p)))
    CSV.write(join([pathstr,"/jcmat","_$(string(excitation_energy)).csv"]), DataFrame(mat, :auto); delim = '\t')
end 