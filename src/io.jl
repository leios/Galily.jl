function write_to_file!(filename, a)
    f = open(filename, "a")
    for i = 1:size(a)[1]
        for j = 1:size(a)[2]
            write(f, string(a[i,j]))
            if j != size(a)[2]
                write(f, '\t')
            end
        end 
        write(f, '\n')
    end 

    write(f, "\n\n")
    close(f)
end
