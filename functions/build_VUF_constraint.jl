function build_vuf_constraint_new(pm_i, vm_var, va_var, vuf_threshold)

    """
    Definition: voltage unbalance is defined as the ratio 
    of the negative sequence voltage component 
    to the positive sequence voltage component
    """

    @variable(pm_i.model, vp_r[1:3])
    @variable(pm_i.model, vp_i[1:3])

    @NLconstraint(pm_i.model, [i = 1:3], vp_r[i] == vm_var[i]*cos(va_var[i])  )
    @NLconstraint(pm_i.model, [i = 1:3], vp_i[i] == vm_var[i]*sin(va_var[i])  )
    

    @variable(pm_i.model, v_pos_r)
    @variable(pm_i.model, v_pos_i)
    @variable(pm_i.model, v_neg_r)
    @variable(pm_i.model, v_neg_i)
    @variable(pm_i.model, v_pos)
    @variable(pm_i.model, v_neg)

    @constraint(pm_i.model, v_pos_r == (1/3)*(vp_r[1] + (-0.5)*vp_r[2] - 0.866*vp_i[2] + (-0.5)*vp_r[3] - (-0.866)*vp_i[3]))
    @constraint(pm_i.model, v_pos_i == (1/3)*(vp_i[1] + (-0.5)*vp_i[2] + 0.866*vp_r[2] + (-0.5)*vp_i[3] + (-0.866)*vp_r[3]))
    @constraint(pm_i.model, v_neg_r == (1/3)*(vp_r[1] + (-0.5)*vp_r[2] - (-0.866)*vp_i[2] + (-0.5)*vp_r[3] - 0.866*vp_i[3]))
    @constraint(pm_i.model, v_neg_i == (1/3)*(vp_i[1] + (-0.5)*vp_i[2] + (-0.866)*vp_r[2] + (-0.5)*vp_i[3] + 0.866*vp_r[3]))
    
    @constraint(pm_i.model, v_pos == v_pos_r^2 + v_pos_i^2)
    @constraint(pm_i.model, v_neg == v_neg_r^2 + v_neg_i^2)


    @constraint(pm_i.model, v_neg <= (vuf_threshold^2)*v_pos)
    @constraint(pm_i.model, v_neg <= (vuf_threshold^2)*v_pos)

end