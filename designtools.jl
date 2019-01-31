
function fourpoint_xy(force, a, l, x, E, I)
    delta = force .*  a ./ (6 .* E .* I) .* (3 .* l .* a .- 3 .* a^2 .- x^2 )
    return delta ./ 1e6 #mm
end

function fourpoint_yx(delta, a, l, x, E, I)
    force = (6 .* delta .* E .* I) ./ (a .* (3 .* l .* a .- 3 .* a^2 .- x^2 ))

    return force .* 1e3 #kN
end

function bearing_res(sec_props, l_bearing,
    fy; E = 200e3, interior = true, phi_bi = 0.80, phi_be = 0.75)
    #input: section, yield strength
    #if interior = true, this is a load point on the INTERIOR of the beam
    #if false, this is a member end bearing support
    h = sec_props.d .- (2 .* sec_props.t_f) #mm, web height
    w = sec_props.t_w #mm, web thickness
    t = sec_props.t_f #mm, flange thickness

    web_slenderness = h ./ w

    if (web_slenderness .>= 1100/sqrt(fy))[1]
        println("Web slenderness limit reached. Bearing stiffeners required.")
    end

    if interior
        println("Interior Load:")
        #interior bearing resistance
        N = l_bearing .+ 10 .* t

        br_1 = phi_bi .* w .* N .* fy / 1e3 #kN, web local plastic buckling
        br_2 = phi_bi .* 1.45 .* w.^2 .* sqrt(fy * E) ./ 1e3 #kN, web overall buck

        println("Web local buckling resistance: " * string(br_1) * " kN")
        println("Web overall buckling resistance: " * string(br_2) * " kN")
        return min(br_1, br_2)

    else
        println("Member end support:")
        #beam end bearing resistance
        N = l_bearing .+ 4 .* t
        br_1 = phi_be .* w .* N .* fy / 1e3 #kN, web local plastic buckling
        br_2 = phi_be .* 0.60 .* w.^2 .* sqrt(fy * E) ./ 1e3 #kN, web overall buck

        println("Web local buckling resistance: " * string(br_1) * " kN")
        println("Web overall buckling resistance: " * string(br_2) * " kN")
        return min(br_1, br_2)
    end
end

function comp_resistance(area, fy, ry, k, L; E = 200e3, phi = 0.90, n = 1.34)
    lambda = k .* L ./ ry .* sqrt.(fy ./ pi^2 ./ E)

    return phi .* area .* fy .* (1 .+ lambda.^(2n)).^(-1 ./ n) ./ 1e3
end

function bearing_stiffener_res(sec_props, t_stiffener, d_stiffener,
        fy; n_stiffeners = 1.0, k = 0.75, E = 200e3, interior = true)

        #check stiffener dimension is within flange width:
        d_max = (sec_props.b_f .- sec_props.t_w) ./ 2

        if (d_stiffener .> d_max)[1]
            println("Bearing stiffeners extend past flange.")
            println("Max. stiffener depth = " * string(d_max) *" mm.")
            return
        end

        #web height
        h_w = sec_props.d .- (2 .* sec_props.t_f)

        #slenderness limit check
        stiff_slenderness = d_stiffener ./ t_stiffener
        if (stiff_slenderness .>= (200 / sqrt(fy)))[1]
            println("Stiffener too slender. Revise.")
            return
        end

        #web area
        if interior
            b_web = 25 .* (sec_props.t_w) .- t_stiffener
            a_web = 25 .* (sec_props.t_w).^2
        else
            b_web = 12 .* (sec_props.t_w) .- t_stiffener
            a_web = 12 .* (sec_props.t_w).^2
        end

        #total effective area
        a_total = a_web .+ (d_stiffener .* t_stiffener) .* 2 .* n_stiffeners

        #moment of inertia
        #stiffener
        d_stiff = 2 .* d_stiffener .+ sec_props.t_w
        b_stiff = t_stiffener

        #web
        d_web = sec_props.t_w
        #b_Web is defined above

        I = n_stiffeners .* (b_stiff .* d_stiff.^3) ./ 12 +
            (b_web .* d_web.^3) ./ 12

        #radius of gyration
        r = sqrt.(I ./ a_total)

        comp = comp_resistance(a_total, fy, r, k, h_w, phi = 1)

        println("Bearing resistance: " * string(comp) * " kN")
        return comp
end

function init_plate_sizer(sec_props, m_p_nominal, fy_p)

    #input the gross plastic moment of a section,
    #returns the required exterior plate thickness to transfer the moment as
    #a force couple.
    #this function ignores the contribution of the interior plates.
    typ_thickness = [6.35, 12.7, 16, 19, 20, 22,
        22.2, 24, 25.4, 27, 28.6, 30, 31.8, 36, 38.1]


    for thickness in typ_thickness
        force = thickness .* sec_props.b_f .* fy_p ./ 1e3 #yield force

        moment = force .* (sec_props.d .+ thickness) ./ 1e3 #kNm

        if (moment .>= m_p_nominal)[1]
            ratio = moment ./ m_p_nominal
            println("Plate thickness of ",
                thickness, "mm provides ",
                ratio .* 100, "% of required gross plastic moment.")
            return thickness
            break
        end
    end
    println("Thickness must exceed 38.1mm.")
    return
end

function init_plate_sizer2(sec_props, d_b, gauge, fy_beam, fy_plate, fu_plate;
        tol = 0, n_rows = 2)

    #input the gross plastic moment of a section,
    #returns the required exterior plate thickness to transfer the moment as
    #a force couple.

    #Typical thicknesses of steel plates
    typ_thickness = [6.35, 12.7, 16, 19, 20, 22,
        22.2, 24, 25.4, 27, 28.6, 30, 31.8, 36, 38.1]

    #determine max gross capacity of beam section
    m_plastic = plastic_moment(sec_props, fy_beam) #kNm

    #force in flange required to be transmitted by plates
    flange_force = moment2force(sec_props, m_plastic) #kN

    #force is split in half between interior and exterior plates
    plate_force = flange_force ./ 2 #kN

    exterior_force = 2 .* plate_force
    interior_force = plate_force

    #design exterior plate
    for thickness in typ_thickness
        ext_gross = thickness .* sec_props.b_f #mm^2
        ext_net = ext_gross .- ((d_b .+ tol) .* thickness) .* n_rows #mm^2

        ext_gross_yield = ext_gross .* fy_plate ./ 1e3
        ext_net_frac = ext_net .* fu_plate ./ 1e3

        ext_critical = min(ext_gross_yield, ext_net_frac)

        if (ext_critical .>= exterior_force)[1]
            ratio = ext_critical ./ exterior_force .* 100
            println("Exterior Plate Thickness of ", thickness, "mm provides ",
                ratio, "% of required force.")
            ext_thickness = thickness
            break
        end
    end


    #design interior plate
    for thickness in typ_thickness
        int_gross = thickness .* (sec_props.b_f ./ 2 .- sec_props.k_1)
        int_net = int_gross .- ((d_b .+ tol) .* thickness) .* (n_rows ./ 2)

        int_gross_yield = int_gross .* fy_plate ./ 1e3
        int_net_frac = int_net .* fu_plate ./ 1e3

        int_critical = min(int_gross_yield, int_net_frac)

        if (int_critical .>= interior_force)[1]
            ratio = int_critical ./ interior_force .* 100
            println("Interior Plate Thickness of ", thickness, "mm provides ",
                ratio, "% of required force.")
            int_thickness = thickness
            break
        end
    end


    return ext_thickness, int_thickness
end
