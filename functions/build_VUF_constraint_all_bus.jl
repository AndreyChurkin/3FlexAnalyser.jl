"""
This function imposes VUF constraints for multiple buses.
The constraints are explicitly included in the optimisation model.
"""

function build_vuf_constraint_allbus(pm_i, vm_var, va_var, vuf_threshold, N_bus)

    """
    Definition: voltage unbalance is defined as the ratio 
    of the negative sequence voltage component 
    to the positive sequence voltage component
    """

    @variable(pm_i.model, vp_r[bus=1:N_bus, i=1:3])
    @variable(pm_i.model, vp_i[bus=1:N_bus, i=1:3])

    @NLconstraint(pm_i.model, [bus=1:N_bus, i=1:3; !(bus in exclude_buses_from_vuf_constraints)], vp_r[bus,i] == vm_var[bus][i]*cos(va_var[bus][i])  )
    @NLconstraint(pm_i.model, [bus=1:N_bus, i=1:3; !(bus in exclude_buses_from_vuf_constraints)], vp_i[bus,i] == vm_var[bus][i]*sin(va_var[bus][i])  )
    

    @variable(pm_i.model, v_pos_r[bus=1:N_bus])
    @variable(pm_i.model, v_pos_i[bus=1:N_bus])
    @variable(pm_i.model, v_neg_r[bus=1:N_bus])
    @variable(pm_i.model, v_neg_i[bus=1:N_bus])
    @variable(pm_i.model, v_pos[bus=1:N_bus])
    @variable(pm_i.model, v_neg[bus=1:N_bus])

    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_pos_r[bus] == (1/3)*(vp_r[bus,1] + (-0.5)*vp_r[bus,2] - 0.866*vp_i[bus,2] + (-0.5)*vp_r[bus,3] - (-0.866)*vp_i[bus,3]))
    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_pos_i[bus] == (1/3)*(vp_i[bus,1] + (-0.5)*vp_i[bus,2] + 0.866*vp_r[bus,2] + (-0.5)*vp_i[bus,3] + (-0.866)*vp_r[bus,3]))
    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_neg_r[bus] == (1/3)*(vp_r[bus,1] + (-0.5)*vp_r[bus,2] - (-0.866)*vp_i[bus,2] + (-0.5)*vp_r[bus,3] - 0.866*vp_i[bus,3]))
    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_neg_i[bus] == (1/3)*(vp_i[bus,1] + (-0.5)*vp_i[bus,2] + (-0.866)*vp_r[bus,2] + (-0.5)*vp_i[bus,3] + 0.866*vp_r[bus,3]))
    
    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_pos[bus] == v_pos_r[bus]^2 + v_pos_i[bus]^2)
    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_neg[bus] == v_neg_r[bus]^2 + v_neg_i[bus]^2)


    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_neg[bus] <= (vuf_threshold^2)*v_pos[bus])
    @constraint(pm_i.model, [bus=1:N_bus; !(bus in exclude_buses_from_vuf_constraints)], v_neg[bus] <= (vuf_threshold^2)*v_pos[bus])

end