*! version 1.0.0  18jul2026
*! nnsreg_estat: clean rejection of unsupported estat calls after nnsreg
*! Author: Eric A. Booth (Senior Researcher, Texas 2036, eric.a.booth@gmail.com)

program define nnsreg_estat
    version 16.0
    if "`e(cmd)'" != "nnsreg" {
        error 301
    }
    di as error "estat is not supported after nnsreg"
    exit 321
end
