#!/bin/sh

genbox << EOF
section.box
EOF
#
reatore2 << EOF
box
incline
EOF
#
rm -f box.rea
rm -f incline.rea
genmap << EOF
incline
0.00001
EOF

