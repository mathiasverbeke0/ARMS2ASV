with open("/home/guest/Internship/data/databases/GEANS_reflib_COI_unaligned3.fas") as data:

    Lines = data.readlines()
    
    with open("/home/guest/Internship/data/databases/GEANS_reflib_COI_unaligned4.fas", 'w') as f:
        
        for line in Lines:
            if line[0] == '>':
                line_split = line.split(';')

                line = ';'.join(line_split[1:len(line_split)])

                f.write(f'>{line}')
                
            
            else:
                if '-' in line:
                    line = line.replace('-', '')
                
                f.write(line)
                    