#!/bin/csh

if ( $#argv != 1 ) then
    echo "Usage: csh script.csh <FRB_name>"
    goto death
endif
set FRB = $1

echo $FRB
cd $FRB/ascii_files/

foreach beam ( 01 02 03 04 05 06 07 08 09 10 11 12 13 )

#    cd $beam
    echo "Beam" $beam

    # Make some outputs in various formats
    pdv -t -j "B 500" ../ar_files/$beam".ar" | awk 'NR!=1{print $2,$3,$4}' > $beam".4ms.cube"
    pdv -t ../ar_files/$beam".ar" | awk 'NR!=1{print $2,$3,$4}' > $beam".raw.cube"
    cd ..

end
cd ../../

death:
exit
