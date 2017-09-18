camera [
eye 0 -17 6
focus 0 0 0
up 0 0 1
angle 60
near-far 0.1 30
]

light [
type point
position -3 4 6
color .8 .8 .8
function 0.5 0 0
]

light [
position 5 -6 -3
color .8 .8 .8
function 0.5 0 0
]

light [
position -5 -5 2.5
color .8 .8 .8
function 0.5 0 0
]

mastersubgraph disk
[
        trans
        [
                scale  5 1 5
                object cylinder
                diffuse 1 1 1
        ]
        ]
        trans
        [
                translate 0 2 0
                scale  2 10 2
                object cylinder
                diffuse 1 1 1
        ]
        ]
]

mastersubgraph root [
trans
        [
                subgraph disk
        ]
trans
        [
                translate 0 -3 0		
                subgraph disk
        ]

trans [
  translate 0 -3 0.75
  scale 5.9 4.3 1.5
  object cube [
    diffuse 0.2 0.1 0.3
  ]
]

trans [
  translate 0 0 0.75
  scale 5.9 4.3 1.5
  object cube [
    diffuse 0.2 0.1 0.3
  ]
]

trans [
	translate 0 8.5 2.9
	scale 6.2 4.5 5.8
	object cube [
		diffuse 0.3 0.2 0.1
	]
]
]
