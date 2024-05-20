def get_name(**kwargs):
    integrator_type = kwargs.get("data_type", "fmcw")
    duration = kwargs.get("T")
    depth = kwargs.get("maxDepth")
    scale = kwargs.get("scale")
    spp = kwargs.get("spp")
    
    output_xml_file_name = "fmcw_correlate_%s_static_duration_%d_depth_%d_scale_%d_spp_%d_%s" % (path_correlation, (int)(T), maxDepth, scale, spp, light_option_name)
    