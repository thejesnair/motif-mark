#!/usr/bin/env python

#conda create -n my_pycairo pycairo
#conda activate my_pycairo

import cairo

surface = cairo.SVGSurface("pycairo_basics.svg", 400, 500)
context = cairo.Context(surface)
# draw a line woo
context.set_line_width(3)
context.move_to(100,58) # (x, y) (0,0)  100 right, 58 down is start
context.line_to(220,58) # line to the right, 120 to the right
context.stroke() # draws a line b/w both points
# rectangle time
context.rectangle(180,115,150,40) #(x0,y0,x1,y1)
context.fill_preserve()
context.set_source_rgb(.5,.4,.6)
context.fill()
surface.finish()




