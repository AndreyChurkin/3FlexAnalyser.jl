"""
This function imposes phase coordination constraints.
That is, only units connected to the phase under consideration are allowed to provide flexibility services.
Other units cannot provide flexibility or manage voltage constraints.
"""

function build_phase_coordination_constraints(pm_i, phase_i)

    for gen_i = 1:length(pm_i.var[:it][:pmd][:nw][0][:pg])
        if gen_i != source_gen_i
            for phase = 1:3
                if phase != phase_i
                    if math["gen"][string(gen_i)]["connections"][1] == phase
                        @constraint(pm_i.model, -10^-5 <= pm_i.var[:it][:pmd][:nw][0][:pg][gen_i][phase] <= 10^-5)
                        @constraint(pm_i.model, -10^-5 <= pm_i.var[:it][:pmd][:nw][0][:qg][gen_i][phase] <= 10^-5)
                    end
                end
            end
        end

    end

end