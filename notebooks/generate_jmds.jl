using Weave
files = readdir(pwd())
filename(path) = splitext(basename(path))[1]

for f in files 
    if occursin(".ipynb", f) && f != ".ipynb_checkpoints"
        convert_doc(f, filename(f) * ".jmd")
    end
end