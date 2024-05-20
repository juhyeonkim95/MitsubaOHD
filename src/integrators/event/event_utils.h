#include <mitsuba/render/scene.h>
#include <mitsuba/core/statistics.h>


using namespace mitsuba;

Float safe_log(Float x, Float threshold=1e-8){
    return x > threshold ? std::log(x) : std::log(threshold);
}

struct PiecewiseLinearSamples
{
    std::vector<Float> times;
    std::vector<Float> values;

    void merge(PiecewiseLinearSamples &other){
        std::vector<Float> new_times;
        std::vector<Float> new_values;
        
        int i = 0, j = 0, k = 0;
        // Traverse both array
        while (i < times.size() && j < other.times.size())
        {
            // Check if current element of first
            // array is smaller than current element
            // of second array. If yes, store first
            // array element and increment first array
            // index. Otherwise do same with second array
            if (times.at(i) < other.times.at(j))
            {
                new_times.push_back(times.at(i));
                new_values.push_back(values.at(i));
                i+= 1;
            } else {
                new_times.push_back(other.times.at(j));
                new_values.push_back(other.values.at(j));
                j+= 1;
            }
        }
    
        // Store remaining elements of first array
        while (i < times.size()){
            new_times.push_back(times.at(i));
            new_values.push_back(values.at(i));
            i+= 1;
        }
    
        // Store remaining elements of second array
        while (j < other.times.size()){
            new_times.push_back(other.times.at(j));
            new_values.push_back(other.values.at(j));
            j+= 1;
        }
        times=new_times;
        values=new_values;
    }

    void add_sample(Float time, Float value){
        times.push_back(time);
        values.push_back(value);
    }

    Float get_value(Float time){
        for(int i=0; i<times.size(); i++){
            if(times.at(i) > time){
                if(i == 0){
                    return values.at(0);
                }
                float dt = (times.at(i) - times.at(i-1));
                // if(dt == 0){
                //     return values.at(i - 1);
                // }
                float r = (time - times.at(i-1)) / dt;
                return values.at(i - 1) * (1 - r) + values.at(i) * r;
            }
        }
        return values.at(values.size() - 1);
    }

    Float get_time(Float value, Float min_time){
        for(int i=0; i<times.size() - 1; i++){
            if(times.at(i) < min_time){
                continue;
            }
            bool inc = values.at(i) < value && value < values.at(i+1);
            bool dec = values.at(i) < value && value < values.at(i+1);
            if(inc || dec){
                float r = (value - values.at(i)) / (values.at(i+1) - values.at(i));
                return times.at(i) * (1 - r) + times.at(i + 1) * r;
            }
        }
        return 10.0;
    }

    void put_to(Float *temp, size_t M){
        int j=0;
        int ne = values.size();
        for(int m=0; m<M; m++){
            float time = m / (float) M;
            if(time < times.at(0)){
                temp[m] = values.at(0);
                continue;
            }
            while(true){
                if(j+1 >=ne){
                    break;
                }
                if (j+1 < ne && time < times.at(j+1)){
                    break;
                }
                j+=1;
            }
            // while((j <= event_times.size() - 2) && (time > event_values.at(j+1))){
            //     j += 1;
            // }
            Float value;
            if(j == times.size() - 1){
                value = values.at(j);
            }
            else{                        
                Float t1 = times.at(j);
                Float t2 = times.at(j + 1);
                Float f1 = values.at(j);
                Float f2 = values.at(j + 1);

                Float r = (time - t1) / (t2 - t1);
                value = r * (f2 - f1) + f1;
            }
            temp[m] = value;
        }
        return;
    }
};

void put_event_data(Float *temp, size_t M, std::vector<float> event_times, std::vector<float> event_values){
    int j=0;
    int ne = event_values.size();
    for(int m=0; m<M; m++){
        float time = m / (float) M;
        while(true){
            if(j+1 >=ne){
                break;
            }
            if (j+1 < ne && time < event_times.at(j+1)){
                break;
            }
            j+=1;
        }
        // while((j <= event_times.size() - 2) && (time > event_values.at(j+1))){
        //     j += 1;
        // }
        Float value;
        if(j == event_times.size() - 1){
            value = event_values.at(j);
        }
        else{                        
            Float t1 = event_times.at(j);
            Float t2 = event_times.at(j + 1);
            Float f1 = event_values.at(j);
            Float f2 = event_values.at(j + 1);

            Float r = (time - t1) / (t2 - t1);
            value = r * (f2 - f1) + f1;
        }
        temp[m] = value;
    }
    return;
}