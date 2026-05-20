using TransverseFol
for i in 4:100
	TransverseFol.runjob(i,2,reg=true)
	TransverseFol.quickview(i,2, obstructions=true, save_html=true, save_png=true,png_width=1920, png_height=1080, png_scale=1, font_size=30)
end