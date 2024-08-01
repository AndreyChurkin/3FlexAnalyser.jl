"""
This function calculates the voltage unbalance factor (OF) a posteriori.
That is, it computes VUF using bus voltages after an OPF problem has been solved.
It can be used for analysing voltage unbalance in case studies and certain OPF solutions.
"""

function VUF_calculation(vm_var, va_var)
        """
        Definition: voltage unbalance is defined as the ratio 
        of the negative sequence voltage component 
        to the positive sequence voltage component
        """

        # va in radians:
        vaa = value.(va_var[1])
        vab = value.(va_var[2])
        vac = value.(va_var[3])

        # get vm value:
        vma = value.(vm_var[1])
        vmb = value.(vm_var[2])
        vmc = value.(vm_var[3])

        # change phasor to complex form:
        vpa = vma*cos(vaa) + (vma*sin(vaa))im
        vpb = vmb*cos(vab) + (vmb*sin(vab))im
        vpc = vmc*cos(vac) + (vmc*sin(vac))im

        a = -0.5 + 0.866im
        a_sqr = -0.5 - 0.866im

        v_pos = (1/3)*(vpa + a*vpb + a_sqr*vpc)
        v_neg = (1/3)*(vpa + a_sqr*vpb + a*vpc)

        v_pos = sqrt(real(v_pos)^2 + imag(v_pos)^2)
        v_neg = sqrt(real(v_neg)^2 + imag(v_neg)^2)

        vuf = v_neg / v_pos

    return vuf
end