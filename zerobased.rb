lineno = 0

File.open("newmatrix.txt", "w") do |output|
  File.open("tempmatrix.mtx", "r") do |input|
    while (line = input.gets)
      # increase line number
      lineno += 1

      # skip header
      if lineno <= 2
	      next
	end
      if lineno == 3
        output.puts line
        next
      end

      # shift numbers
      cols = line.split(' ')
      # puts cols.inspect
      output.puts("#{cols[0].to_i-1} #{cols[1].to_i-1} #{cols[2]}")
    end
  end
end
