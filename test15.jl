dict1=load("loc4.jld2")

uber1=dict1["uber"]
cber1=dict1["cber"]

uber_com = dict1["uber_com"]
cber_com = dict1["cber_com"]

plot(uber1)
plot!(cber1)
plot!(uber_com)
plot!(cber_com)
(plot!(title=("̄ber: (12 km range)")))
(plot!(ylabel = ("̄ber")))
(plot!(xlabel = ("packet no.")))
(plot!(xtickfont=font(10)))
(plot!(ytickfont=font(10)))
(plot!(labelfontsize=15))
plot!(ylims=(0.0,0.5))
display(plot!(legend=true))


# dict1=load("loc5.jld2")

# uber1=dict1["uber"]
# cber1=dict1["cber"]

# uber_com = dict1["uber_com"]
# cber_com = dict1["cber_com"]

# plot(uber1)
# plot!(cber1)
# plot!(uber_com)
# plot!(cber_com)
# (plot!(title=("̄ber: (15 km range)")))
# (plot!(ylabel = ("̄ber")))
# (plot!(xlabel = ("packet no.")))
# (plot!(xtickfont=font(10)))
# (plot!(ytickfont=font(10)))
# (plot!(labelfontsize=15))
# plot!(ylims=(0.0,0.4))
# display(plot!(legend=true))

# dict1=load("loc6.jld2")

# uber1=dict1["uber"]
# cber1=dict1["cber"]

# uber_com = dict1["uber_com"]
# cber_com = dict1["cber_com"]

# plot(uber1)
# plot!(cber1)
# plot!(uber_com)
# plot!(cber_com)
# (plot!(title=("̄ber: (17 km range)")))
# (plot!(ylabel = ("̄ber")))
# (plot!(xlabel = ("packet no.")))
# (plot!(xtickfont=font(10)))
# (plot!(ytickfont=font(10)))
# (plot!(labelfontsize=15))
# plot!(ylims=(0.0,0.4))
# display(plot!(legend=true))