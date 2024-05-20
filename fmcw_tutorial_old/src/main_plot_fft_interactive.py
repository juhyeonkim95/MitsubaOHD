import os
import numpy as np
import matplotlib.pyplot as plt
from utils.image_utils import *
from utils.common_path import *
import itertools
from scipy.fft import fft, ifft, fftfreq, fftshift
from utils.exp_configs_utils import *
from skimage.transform import resize

keep_position = False
use_specific = False



def visualize_transient(**kwargs):
    # all of the different options
    scene_name_org = kwargs.get("scene_name")
    scene_names = kwargs.get("scene_names")
    data_options = kwargs.get("data_options")

    use_crop = kwargs.get("use_crop", False)
    use_normalize = kwargs.get("use_normalize", True)

    if scene_name_org is None:
        cropOffsetX, cropOffsetY = get_crop_pos(scene_names[0])
    else:
        cropOffsetX, cropOffsetY = get_crop_pos(scene_name_org)
    cropWidth = 1
    cropHeight = 1

    file_names = kwargs.get("file_names")
    scene_scale = kwargs.get("scene_scale", 1)
    
    xmin, xmax = kwargs.get("xrange")
    scale = kwargs.get("scale", 8)
    scales = kwargs.get("scales")
    show_autocorrelation = kwargs.get("show_autocorrelation", False)

    # FMCW related
    T = kwargs.get("T", 10) * 1e-6
    B = kwargs.get("B", 1) * 1e9
    wavelength = kwargs.get("wavelength", 1550) * 1e-9
    c  = 3e8
    f_c = c / wavelength

    image_output_folder = kwargs.get("image_output_folder")
    
    images = {}
    first_image = None
    bin_size = 2048
    M = 4096
    configs = []
    data_options_new = []

    if file_names is not None:
        for i, file_name in enumerate(file_names):
            data_option = data_options[i]

            if scene_names is not None:
                scene_name_i = scene_names[i]
            else:
                scene_name_i = scene_name_org
            if scene_scale != 1 and data_option != "transient" and data_option != "ground_truth":
                scene_name = "%s_scale_%d" % (scene_name_i, scene_scale)
            else:
                scene_name = scene_name_i

            if data_option == "ground_truth":
                exp_output_folder = os.path.join(image_output_folder, "gt_depth", scene_name)
                output_file_name = os.path.join(exp_output_folder, "depth")
                depth_map_all = np.load(output_file_name+".npy")
                image = depth_map_all[:, :, 0] * scene_scale
                image = resize(image, (image.shape[0] // scale, image.shape[1] // scale))
                #print(image.shape)
                #exit(0)
            else:
                image = np.load(os.path.join(image_output_folder, data_option, scene_name, file_name+".npy"))
                if scales is not None and scales[i] != 8:
                    image = resize(image, (image.shape[0] * scales[i] // 8, image.shape[1] * scales[i] // 8))
                    print("RESHAPED", image.shape)

            if "fmcw_dynamic" in data_option:
                
                image_key1 = "%s_%s_up" % (scene_name, file_name)
                images[image_key1] = image[:,:,0:M]
                image_key2 = "%s_%s_down" % (scene_name, file_name)
                images[image_key2] = image[:,:,M:2*M]
                
                configs.append(image_key1)
                configs.append(image_key2)
                data_options_new.append(data_option)
                data_options_new.append(data_option)
                
                if first_image is None or data_option == "transient":
                    # first_image = image

                    yf1 = fft(image[:,:,0:M], axis=-1)
                    magnitude1 = np.abs(yf1)
                    magnitude_max_idx1 = np.argmax(magnitude1, axis=-1)

                    yf2 = fft(image[:,:,M:2*M], axis=-1)
                    magnitude2 = np.abs(yf2)
                    magnitude_max_idx2 = np.argmax(magnitude2, axis=-1)
                    
                    xf = fftfreq(M, 1) * M
                    fp1 = np.abs(xf[magnitude_max_idx1])
                    fp2 = np.abs(xf[magnitude_max_idx2])

                    first_image = ((fp1 + fp2) * 0.5)
                    first_image = first_image[:,:,None]
            else:
                image_key = "%s_%s" % (scene_name, file_name)
                images[image_key] = image
                if "transient" in data_option:
                    bin_size = image.shape[2]
                configs.append(image_key)
                if first_image is None or data_option == "transient":
                    first_image = image
                data_options_new.append(data_option)
                
                print(image.shape)
                print(image.sum())
                if "transient" in data_option and show_autocorrelation:
                    images[image_key+"_auto"] = image
                    data_options_new.append(data_option+"_auto")
                    configs.append(image_key+"_auto")

                    print(image.shape)
                    print(image.sum())



    data_options = data_options_new
    fig = plt.figure(figsize=kwargs.get("figsize", [8, 4]))

    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    if kwargs.get("show_label", False):
        ax2.set_xlabel("distance")
        ax2.set_ylabel("magnitude")

    data_plots = {}
    styles = kwargs.get("styles")
    
    idx = 0

    for i, config in enumerate(configs):
        style = styles[idx % len(styles)]
        if kwargs.get("labels") is not None:
            label = kwargs.get("labels")[i]
        else:
            label = str(config)
        if "ground_truth" in config:
            data_plot = ax2.axvline(label="gt", color='black', linestyle='--', linewidth=2, alpha=1.0)
        else:
            data_plot, = ax2.plot([], [], label=label, **style)
        data_plots[config] = data_plot
        idx += 1



    depth_image = np.mean(first_image, axis=-1)
    depth_image = to_ldr_image(depth_image)
    
    
    # noises = {}
    # for config, value in images.items():
    #    noise = np.random.random(value.shape[2])
    #    noises[config] = noise

    def get_value(x, y, config):
        
        if file_names is None:
            _, _, _, _, data_option = config
        else:
            index = configs.index(config)
            data_option = data_options[index]
        
        if data_option == "ground_truth":
            depth = images[config][y, x] * 2
            return depth, None
        
        if data_option == "transient":
            data = images[config][y, x, :]
            tmin, tmax = get_transient_min_max_range(config)

            xs = np.linspace(tmin, tmax, bin_size)
            xs *= scene_scale
            return data, xs
        elif data_option == "transient_auto":
            data = images[config][y, x, :]
            tmin, tmax = get_transient_min_max_range(config)

            xs = np.linspace(0, tmax-tmin, bin_size)
            xs *= scene_scale
            auto_corr = np.correlate(data, data, mode='full')
            auto_corr = auto_corr[auto_corr.shape[0] // 2:]
            return auto_corr, xs

        elif "fmcw" in data_option:
            if use_crop and "fmcw" in data_option:
                if x >= cropOffsetX and x < cropOffsetX + cropWidth and y >= cropOffsetY and y < cropOffsetY + cropHeight:
                    data = images[config][y-cropOffsetY, x-cropOffsetX, :]
                else:
                    data = np.ones(images[config].shape[2], )
            else:
                data = images[config][y, x, :]

            # print(images)
            # print(np.sum(data))
            # data = np.nan_to_num(data)
            # print(np.sum(data))
            data = np.array(data)
            # data += noises[config] * 1e-4

            # if data_option == "fmcw_dynamic":
            #     yf1 = fft(data[0:M])
            #     yf1 = fftshift(yf1)
            #     yf1 = np.abs(yf1)
                
            #     #print(yf1.shape, "AA")
            #     #exit(0)

            #     yf2 = fft(data[M:2*M])
            #     yf2 = fftshift(yf2)
            #     yf2 = np.abs(yf2)
                
            #     magnitude = (yf1+yf2)*0.5

            # else:
                # print(data.shape, data_option)

            yf = fft(data)
            yf = fftshift(yf)

            magnitude = np.abs(yf)
            phase = np.angle(yf)

            xf = fftfreq(M, T)
            xf = fftshift(xf) * M
            xf = xf * c * T / (2 * B) * 2
            return magnitude, xf


    def update(x, y):
        for config in configs:
            ydata, xdata = get_value(x, y, config)
            if xdata is None:
                data_plots[config].set_xdata(ydata)
                continue
            
            # normalize
            
            if use_normalize and np.max(ydata) > 0:
                ydata = ydata / np.max(ydata)

            if kwargs.get("use_log_scale", False):
                # ydata += np.power(10, -1.3)
                ydata += np.power(10, -2.0)
                ydata = 20 * np.log10(ydata)
        
            data_plots[config].set_ydata(ydata)
            data_plots[config].set_xdata(xdata)
            
            ax2.set_xlim(xmin = xmin * scene_scale, xmax = xmax * scene_scale)

        if kwargs.get("use_log_scale", False):
            # ax2.set_ylim(ymin = -20.0, ymax = 0.0)
            ax2.set_ylim(ymin = -40.0, ymax = 0.0)
        else:
            ax2.set_ylim(ymin = 0.0, ymax = 1.0)

        if not use_normalize:
            ax2.set_ylim(ymin = np.min(ydata), ymax = np.max(ydata))
        
        if kwargs.get("show_legend", False):
            ax2.legend()
        fig.canvas.draw()

    def mouse_move(event):
        x, y = event.xdata, event.ydata
        if x is None or y is None or keep_position:
            return
        if not use_specific:
            update(int(x), int(y))
        else:
            update(cropOffsetX, cropOffsetY)

    def on_press(event):
        global keep_position, use_specific
        if event.key == 'j':
            keep_position = not keep_position
        elif event.key == 'h':
            use_specific = not use_specific


    plt.connect('key_press_event', on_press)
    plt.connect('motion_notify_event', mouse_move)

    ax1.imshow(depth_image)
    ax1.set_title("depth map")
    ax1.axis('off')

    
    if kwargs.get("save_file", False):
        import matplotlib.transforms as mtransforms
        bbox_inches=mtransforms.Bbox(
            # This is in "figure fraction" for the bottom half
            # input in [[xmin, ymin], [xmax, ymax]]
            [[0.5, 0], [1.0, 1.0]]
        ).transformed(
            # this take data from figure fraction -> inches
            #    transFigrue goes from figure fraction -> pixels
            #    dpi_scale_trans goes from inches -> pixels
            (fig.transFigure - fig.dpi_scale_trans)
        )
        update(cropOffsetX, cropOffsetY)
        output_folder = os.path.join(kwargs.get("plot_output_path"), scene_name)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        plt.tight_layout()
        file_name = kwargs.get("plot_name", "plot")
        #plt.axis('off')
        #ax2.set_xticks([])
        #ax2.set_yticks([])
        #ax2.get_xaxis().set_visible(False)
        #ax2.get_yaxis().set_visible(False)
        extent = ax2.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        plt.savefig(os.path.join(output_folder, "%s.png"%(file_name)), bbox_inches=extent.expanded(1.1, 1.2))
        ax2.set_xticklabels([])
        ax2.set_yticklabels([])
        plt.savefig(os.path.join(output_folder, "%s.svg"%(file_name)), dpi=600, bbox_inches=extent.expanded(1.1, 1.2))
        plt.close('all')
    else:
        plt.show()


def plot_fmcw_dynamic_result():
    scene_name = "cornell-box"

    data_options = [
        "fmcw_dynamic"
    ]
    
    file_names = [
        "fmcw_dynamic_duration_10_depth_2_scale_8_spp_1_collimated_amplitude"
    ]

    styles = [
        {"linewidth": 2, "alpha":1, "color": "C0"},
        {"linewidth": 2, "alpha":1, "color": "C1", "linestyle":"--"},
    ]

    labels = [
        "up-chirp",
        "down-chirp",
    ]

    configs = {
        "scene_name": scene_name,
        "data_options": data_options,
        "image_output_folder": os.path.join(PROJECT_PATH, "results"),
        "file_names":file_names,
        "styles": styles,
        "scene_scale":10,
        "xrange": (0, 35),
        "labels": labels,
        "show_legend": True
    }
    visualize_transient(
        **configs
    )

if __name__ == "__main__":
    plot_fmcw_dynamic_result()
    exit(0)