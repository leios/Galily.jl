function project_to_3d(a; l=1)
    output = zeros(size(a)[1],3)
    for i = 1:size(a)[1]
        temp_w = 1/(l-a[i,4])
        proj_array = [temp_w 0 0 0;
                      0 temp_w 0 0;
                      0 0 temp_w 0]
        output[i,:] .= proj_array*a[i,:]
    end

    return output
end

function write_to_file!(filename, a; project = false, l=1)
    f = open(filename, "a")

    if size(a)[2] == 4 && project
        a = project_to_3d(a, l=l)
    elseif size(a)[2] > 4
        println("We do not support arrays greater than 4D for projection!")
    end

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

