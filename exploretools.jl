using Plots
using DataFrames
using CSV
using LaTeXStrings
pyplot()

W_data = CSV.read("W.csv", allowmissing = :none)

function name2data(section_name)
    #input a section name as a string (eg. "W460X19"),
    #returns all section properties of that section
    #as a database to then be called
    out_frame = W_data[W_data[1:end, :name] .== section_name, :]
    out_frame = out_frame[3:end]
    return out_frame
end

function plastic_moment(sec_props, yield_strength; phi = 1.0)
    #input is one row dataframe of a given section + material yield strength
    #returns the theoretical plastic moment capacity (Z x Fy)
    return float(sec_props.Z_x  * yield_strength ./ 1e3)
end

function s16_14_1(sec_props, d_b, fy; n_rows = 2, tolerance = 0)
    #input section dataframe + bolt diameter in mm, tolerance is default 0
    gross_flange_area = sec_props.b_f .* sec_props.t_f

    area_loss = sec_props.t_f .* (d_b + tolerance) .* n_rows

    net_area = gross_flange_area .- area_loss
    ratio = net_area ./ gross_flange_area

    if (ratio .< 0.85)[1]
        println("85% rule triggered. Net properties should be used.")
        net_zx = net_plastic_modulus(sec_props, d_b, n_rows)
        net_plastic_moment = net_zx .* fy ./ 1e3 #kNm
        println("Net Zx = ", net_zx, "x10^3 mm^3 ; Mp_net = ",
            net_plastic_moment, "kNm.")
        return true, ratio
    else
        println("85% rule not triggered. Gross section properties OK.")
        return false, ratio
    end
end

function edge_distance(sec_props, d_b, gauge)
    #returns edge distance of given section + connection gauge
    edge = (sec_props.b_f .- gauge) ./ 2
    return float(edge)
end

function geo_limits(sec_props, d_b, n_bolts)
    #where n_bolts is IN THE DIRECTION OF LOADING
    #returns a 1x3 ARRAY of minimum pitch, edge and end distances in mm

    pitch_min = 2.7 * d_b #mm

    #minimum edge distances corresponding to bolt diameter (assuming imperial)
    bolt_sizes = 25.4 .* [5/8 3/4 7/8 1 1.125 1.25] #mm
    edge_limits = [28 32 38 44 51 57] #mm (Use conservative sheared edge values)
    idx = findall(x -> x == d_b, bolt_sizes)
    edge_min = edge_limits[idx]

    #minimum end distance
    if n_bolts > 2
        end_min = edge_min
    else
        end_min = 1.5 * d_b
    end

    # println("Minimum values for [pitch edge end] distances (mm)")
    return [pitch_min edge_min end_min]
end

function geo_comparison(design, minimums)
    #input is: [PITCH EDGE END]
    #compares design geometry to s16 limits

    comp = design .> minimums #outputs array of booleans

    pitch_check = comp[1]
    edge_check = comp[2]
    end_check = comp[3]

    if all(value -> value == true, comp)
        println("Bolt hole geometry limits satisfied. Design OK.")
        return true
    else
        println("Bolt hole geometry limits not satisfied. Revise.")
        println("Minimum [pitch edge end] = ", minimums, " ; Design = ", design)
        return false
    end
end

function net_flange_area(sec_props, d_b; n_rows = 2, tolerance = 0)
    gross = sec_props.b_f .* sec_props.t_f #mm^2
    loss = sec_props.t_f .* (d_b + tolerance) .* n_rows #mm^2

    return float(gross .- loss) #mm^2
end

function flange_fracture_webyield(sec_props,
        end_dist, pitch, n_bolts, net_flange_area, fy, fu; phi = 1.0, phi_u = 1.0)

    fracture = phi_u .* net_flange_area .* fu / 1e3
    yield = phi .* (end_dist .+ pitch .* (n_bolts .- 1)) .* sec_props.t_w .* fy ./ 1e3

    return fracture .+ yield
end

function flange_fracture2(sec_props,
        end_dist, pitch, n_bolts, net_flange_area, fy, fu; phi = 1.0, phi_u = 1.0)

    fracture = phi_u .* net_flange_area .* fu / 1e3

    return fracture
end

function s16_block_strength(sec_props, d_b, fy, fu, end_dist,
    pitch, gauge, n_bolts; disp = true, phi_u = 1.0, ut = 1.0, tolerance = 0)

    #d_b = bolt diameter
    #fy, fu = material strengths
    #end_dist = dist from end center of hole to free end of section
    #edge_dist = dist from center of hole to section edge
    #pitch = distance between bolts parallel to load
    #gauge - distance between bolts perp. to load (across the web)
    #n_bolts = number of bolts in one line

    #Edge fracture
    net_fracture1 = (sec_props.b_f .- gauge .- (d_b + tolerance)) .* sec_props.t_f
    gross_shear1 = 2 * (end_dist + pitch * (n_bolts - 1)) .* sec_props.t_f

     #Center fracture
    net_fracture2 = (gauge .- (d_b + tolerance)) .* sec_props.t_f #mm^2
    gross_shear2 = 2 * (end_dist + pitch * (n_bolts - 1)) .* sec_props.t_f
        .+ (end_dist + pitch * (n_bolts - 1)) .* sec_props.t_w

    #gross shear strength
    if fy > 460
        shear_strength = fy
    else
        shear_strength = 0.5 * (fy + fu)
    end

    tension1 = (ut * fu) .* net_fracture1 ./ 1e3
    tension2 = (ut * fu) .* net_fracture2 ./ 1e3

    shear1 = gross_shear1 .* shear_strength .* 0.6 ./ 1e3
    shear2 = gross_shear2 .* shear_strength .* 0.6 ./ 1e3

    block1 = tension1 .+ shear1
    block2 = tension2 .+ shear2

    if (block1 .< block2)[1]
        if disp
            println("Edge block shear governs.")
        end
        return (phi_u .* block1)[1] #kN
    else
        if disp
            println("Center block shear governs.")
        end
        return (phi_u .* block2)[1] #k
    end
end

function force2moment(sec_props, flange_force)
    #flange force should be in kN
    return float((sec_props.d - sec_props.t_f)[1] .* flange_force[1] ./1e3) #kNm
end

function moment2force(sec_props, moment)
    #moment in kNm
    return float((moment .* 1e3) ./ (sec_props.d .- sec_props.t_f)) #kN
end

function block_shear_components(sec_props, d_b, fy, fu, end_dist, edge_dist,
    pitch, gauge, n_bolts ; phi_u = 1.0, ut = 1.0, tolerance = 0)

    #returns the shear and tension components of block shear capacity
    #as simple forces

    #Edge fracture
    net_fracture1 = (sec_props.b_f .- gauge .- (d_b + tolerance)) .* sec_props.t_f
    gross_shear1 = 2 * (end_dist + pitch * (n_bolts - 1)) .* sec_props.t_f

     #Center fracture
    net_fracture2 = (gauge .- (d_b + tolerance)) .* sec_props.t_f #mm^2
    gross_shear2 = 2 * (end_dist + pitch * (n_bolts - 1)) .* sec_props.t_f
        .+ (end_dist + pitch * (n_bolts - 1)) .* sec_props.t_w

    #gross shear strength
    if fy > 460
        shear_strength = fy
    else
        shear_strength = 0.5 * (fy + fu)
    end

    tension1 = (ut * fu) .* net_fracture1 ./ 1e3
    tension2 = (ut * fu) .* net_fracture2 ./ 1e3

    shear1 = gross_shear1 .* shear_strength .* 0.6 ./ 1e3
    shear2 = gross_shear2 .* shear_strength .* 0.6 ./ 1e3

    block1 = tension1 .+ shear1
    block2 = tension2 .+ shear2

    if (block1 .< block2)[1]
        println("Edge block shear governs.")
        tension = tension1 ./ 2 #kN
        shear = shear1 ./ 2 #kN

        println("Tension Capacity = ", tension, "kN (per block)")
        println("Shear Capacity = ", shear, "kN (per block)")

        return tension, shear
    else
        println("Center block shear governs.")
        tension = tension2 ./ 2
        shear = shear2 ./ 2

        println("Tension Capacity = ", tension, "kN")
        println("Shear Capacity = ", shear, "kN (per block)")

        return tension, shear
    end
end


function net_moment_inertia(sec_props, d_b, n_bolts; tol = 0)
    I_initial = sec_props.I_x .* 1e6 #mm^4

    I_bolt = (d_b + tol) .* sec_props.t_f.^3 ./ 12 #mm^4
    A_bolt = (d_b + tol) .* sec_props.t_f #mm^2

    bolt_to_centroid = (sec_props.d .- sec_props.t_f) ./ 2

    I_reduction = n_bolts .* (I_bolt .+ (A_bolt .* bolt_to_centroid.^2))

    I_final = I_initial .- I_reduction

    return I_final ./ 1e6
end

function net_plastic_modulus(sec_props, d_b, n_bolts; tol = 0)

    Z_nominal = sec_props.Z_x

    #base on simplified plastic modulus equation (ignoring contribution of K area)

    b = sec_props.b_f .- (d_b .+ tol) .* n_bolts
    tf = sec_props.t_f
    d = sec_props.d
    tw = sec_props.t_w

    Z_final = b .* tf .* (d .- tf) + 0.25 .* tw .* (d .- tf).^2 #mm^3


    return Z_final ./ 1e3
end

function section_overview(section_name, d_b, end_dist, fy, fu,
    pitch, gauge, n_bolts; phi = 1.0, phi_u = 1.0, tol = 0, disp_plots = true)
    println("For " * section_name * ": ")
    section = name2data(section_name)
    #output theoretical full section moment resistance
    mp_gross = plastic_moment(section, fy)
    println("Gross plastic moment capacity = " * string(mp_gross) * "kNm")

    #check on 85% rule
    redux, ratio = s16_14_1(section, d_b, fy)
    #if redux = true, 85% rule is initiated (IE net section < 0.85gross)

    #check geo limits
    edge = edge_distance(section, d_b, gauge)
    #array of relevant connection geometry
    connex_geo = [pitch edge end_dist]
    #array of minimum distances required by s16 standard
    connex_limits = geo_limits(section, d_b, 2)
    #evaluation of design and required connection geometry
    geo_pass = geo_comparison(connex_geo, connex_limits)

    #flange properties
    net_area = net_flange_area(section, d_b)
    gross_area = section.b_f .* section.t_f

    #flange resistances
    net_fracture = flange_fracture_webyield(section, end_dist, pitch,
        n_bolts, net_area, fy, fu)
    net_fracture2 = flange_fracture2(section, end_dist, pitch,
        n_bolts, net_area, fy, fu)

    gross_yield = gross_area .* fy ./ 1e3
    blockshear_16 = s16_block_strength(section, d_b, fy,
        fu, end_dist, pitch, gauge, n_bolts)

    flange_resistances = [net_fracture gross_yield blockshear_16]

    #moment resistances
    net_fracture_moment = force2moment(section, net_fracture)
    println("Net Fracture Moment (+ web yield) = "
        * string(net_fracture_moment) * "kNm")
    net_fracture_moment2 = force2moment(section, net_fracture2)
    println("Net Fracture Moment (flange only) = "
        * string(net_fracture_moment2) * "knm")
    gross_yield_moment = force2moment(section, gross_yield)
    println("Gross Yield Moment = " * string(gross_yield_moment) * "kNm")
    blockshear_moment = force2moment(section, blockshear_16)
    println("Block Shear Moment = " * string(blockshear_moment) * "kNm")
    #output results
    connex_length = end_dist + pitch * (n_bolts - 1)

    ######
    output_dataframe = DataFrame(Section = section_name,
        BoltDiameter = d_b,
        HoleTolerance = tol,
        EndDistance = end_dist,
        Pitch = pitch,
        Gauge = gauge,
        ConnectionLength = connex_length,
        F_y = fy,
        F_u = fu,
        ϕ = phi,
        ϕ_u = phi_u,
        s16_14_1_trigger = redux,
        NetGrossRatio = ratio,
        GrossPlasticMoment = mp_gross,
        NetFractureMoment1 = net_fracture_moment,
        NetFractureMoment2 = net_fracture_moment2,
        GrossYieldMoment = gross_yield_moment,
        BlockShearMoment = blockshear_moment)

    #Bar plots for moments

    #concatenate the different section moment capacities
    moment_recap = vcat(output_dataframe.GrossPlasticMoment,
        output_dataframe.NetFractureMoment1,
        output_dataframe.NetFractureMoment2,
        output_dataframe.GrossYieldMoment,
        output_dataframe.BlockShearMoment)

    #labels for each moment type
    moment_names = ["Gross Plastic",
    "Net Fracture (Flange+Web Yield)",
    "Net Fracture (Flange Only)",
    "Gross Yield",
    "Block Shear"]

    #display plot
    if disp_plots
        display(bar(moment_names,moment_recap,
            xrotation = 15,
            ylabel = "Moment (kNm)",
            title = "Section Moment Capacity Summary: " * section_name,
            legend = false))
    end

    #blank line for style
    println()

    #returns single row dataframe of chosen section + relevant strengths
    return output_dataframe

    println()
end
