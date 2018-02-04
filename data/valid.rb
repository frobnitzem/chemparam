# Ruby script to validate mol file format
# suitable for sanitizing web interfaces.
# David M. Rogers, Feb. 3, 2018

def valid_mol?(s)
    lines = s.split "\n"
    if lines.any? { |l| l.length > 80 }
        return false, "Line exceeds 80 characters."
    end

    # Header block:
    title_excl = [ /\$MDL/, /\$\$\$\$/, /\$RXN/, /\$RDFILE/ ]
    if title_excl.any? { |m| m =~ lines[0] }
        return false, "Invalid title."
    end
    lines[1] = lines[1][0..22]
    t2 = /^[ a-zA-Z0-9]{10}[[:digit:]]{10}[23]D$/ =~ lines[1]
    if t2.nil?
        return false, "Invalid header."
    end
    lines[2] = "" # sanitize comments

    # Counts line:
    istr  = /^ *[0-9]+$/
    #t3 = /^[ 0-9]{12}  [01][ 0-9]{18} V2000$/ =~ lines[3]
    unless [ lines[3][0..2] =~ istr,
             lines[3][3..5] =~ istr,
             lines[3][12..14] =~ /  [01]/,
             lines[3][33..38] == ' V2000' ].all?
        return false, "Invalid counts line."
    end
    lines[3][6..11] = "  0  0"
    lines[3][15..32] = "  0  0  0  0  0999" # sanitize obsolete fields

    atoms = Integer(lines[3][0..2])
    bonds = Integer(lines[3][3..5])
    #p atoms, bonds
    unless lines.count == 5+atoms+bonds
        return false, "Invalid number of lines."
    end

    # Atoms:
    def isatom?(a)
        tdot4 = /^ *-?[0-9]+\.[0-9]{4}$/
        if a[44] == '0'
            a[44] = '1'
        end
        a[48..-1] = ''
        [ a[0..9] =~ tdot4, a[10..19] =~ tdot4, a[20..29] =~ tdot4,
          a[30..33] =~ /^ [a-zA-Z] *$/,
          a[34..47] =~ /^[- ][0-4]  [0-7]  [0-3]  [1-5]  [0-1]$/
        ].all?
    end
    unless lines[4..(3+atoms)].all? { |a| isatom?(a) }
        return false, "Invalid atom line."
    end

    def isbond?(b)
        istr  = /^ *[0-9]+$/
        if b.length == 18
            b = b + '  0'
        end
        [ b[0..2] =~ istr, b[3..5] =~ istr,
          b[6..20] =~ /^  [1-3]  [01346]  0  [0-2]  0$/,
        ].all?
    end
    unless lines[(4+atoms)..(3+atoms+bonds)].all? { |b| isbond?(b) }
        return false, "Invalid bond line."
    end

    if (/^M  END$/ =~ lines[lines.count-1]).nil?
        return false, "Unterminated mol (requires M  END line)"
    end

    lines = lines[0..(3+atoms+bonds)]
    [true, lines.join("\n") + "\nM  END\n"]
end

f = open("test.mol")
r = f.read
f.close()

ok, msg = valid_mol?( r )
if !ok
    p "Error!"
    p msg
else
    print msg
end

