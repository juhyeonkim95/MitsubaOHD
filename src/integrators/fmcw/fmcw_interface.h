using namespace mitsuba;

class FMCWInterface {
public:
    FMCWInterface() {};
    FMCWInterface(const Properties &props) {
        // m_f_c = props.getFloat("f_c", 77); 
        if(props.hasProperty("wavelength")){
            Float wavelength = props.getFloat("wavelength", 1550); // nm (1e-9)
            m_f_c = 3e8 / wavelength;
        } else{
            m_f_c = props.getFloat("f_c", 77); // GHz (1e9)
        }
        m_B = props.getFloat("B", 1);  // GHz (1e9)
        m_T = props.getFloat("T", 10);  // micro second (1e-6)
        m_time = props.getFloat("time", 0); // micro second (1e-6)
        m_R_min = 2 * props.getFloat("R_min", 0); // m
        m_M = props.getInteger("M", 4096);  // number of sampled time
        m_use_collimated_light = props.getBoolean("use_collimated", false);
        m_fov_error = props.getFloat("fov_error", 0.5);
        m_use_amplitude = props.getBoolean("use_amplitude", false);
        m_sqrt_pdf = props.getBoolean("pdf_sqrt", false);
        bool invert = props.getBoolean("invert", false);
        if(invert){
            m_f_c += m_B;
            m_B = -m_B;
        }


        printf("[FMCWInterface INFO - static]\n");
        printf("\tT (us): %f\n", m_T);
        printf("\tB (GHz): %f\n", m_B);
        printf("\tf_c (GHz): %f\n", m_f_c);
        printf("\ffov_error : %f\n", m_fov_error);
        if(m_use_collimated_light){
                printf("\tm_use_collimated_light: True\n");
        } else { printf("\tm_use_collimated_light: False\n");}

        if(m_use_amplitude){
                printf("\tm_use_amplitude: True\n");
        } else { printf("\tm_use_amplitude: False\n");}
        if(m_sqrt_pdf){
                printf("\tm_sqrt_pdf: True\n");
        } else { printf("\tm_sqrt_pdf: False\n");}
    };

    Float get_fmcw_weight(Float t, Float path_length) const {
        // R / c
        // R / 3e8 * 1e6
        Float t_d = (path_length - m_R_min) / 3e2; 

        // 2 * pi * f_c * t_d
        // 2 * pi * f_c * 1e9 * t_d * 1e-6
        Float A = 2 * M_PI * m_f_c * t_d * 1e3;

        // pi * B * t_d * t_d / T
        // pi * B * 1e9 * t_d * 1e-6 * t_d * 1e-6 / (T * 1e-6)
        Float B = M_PI * m_B * t_d * t_d / m_T * 1e3;

        // 2 * pi * B * t_d * t / T
        // 2 * pi * B * 1e9 * t_d * 1e-6 * t * 1e-6 / (T * 1e-6)
        Float C = 2 * M_PI * m_B * t_d * t / m_T * 1e3;

        return (std::cos(A + B + C));
    }

    Float get_fmcw_weight_invert(Float t, Float path_length) const {
        Float n_f_c = m_f_c + m_B;
        Float n_B = -m_B;

        // R / c
        // R / 3e8 * 1e6
        Float t_d = (path_length - m_R_min) / 3e2; 

        // 2 * pi * f_c * t_d
        // 2 * pi * f_c * 1e9 * t_d * 1e-6
        Float A = 2 * M_PI * n_f_c * t_d * 1e3;

        // pi * B * t_d * t_d / T
        // pi * B * 1e9 * t_d * 1e-6 * t_d * 1e-6 / (T * 1e-6)
        Float B = M_PI * n_B * t_d * t_d / m_T * 1e3;

        // 2 * pi * B * t_d * t / T
        // 2 * pi * B * 1e9 * t_d * 1e-6 * t * 1e-6 / (T * 1e-6)
        Float C = 2 * M_PI * n_B * t_d * t / m_T * 1e3;

        return (std::cos(A + B + C));
    }

protected:
    Float m_time;
    Float m_B;
    Float m_T;
    Float m_f_c;
    Float m_R_min;
    Float m_fov_error;
    size_t m_M;
    bool m_use_collimated_light;
    bool m_use_amplitude;
    bool m_sqrt_pdf;
};