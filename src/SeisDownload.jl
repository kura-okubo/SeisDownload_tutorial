__precompile__()
module SeisDownload

include("utils.jl")
using .Utils, SeisIO, Dates, Printf, JLD2, FileIO,  Distributed

export ParallelSeisrequest, seisdownload_NOISE

#------------------------------------------------------------------#
#For the time being, we need remove_response function from obspy
#This will be replaced by SeisIO modules in the near future.
#Please activate obspy enviroment before launching Julia.
include("remove_response_obspy.jl")
using .Remove_response_obspy
#------------------------------------------------------------------#


"""
    ParallelSeisrequest(NP::Int, InputDict::Dict)

    Request seismic data with Multiple cores.
# Arguments
- `NP`           : number of processors
- `InputDict`    : dictionary which contains request information
"""
function ParallelSeisrequest(NP::Int, InputDict::Dict)

    Utils.initlogo()

    #stationlist
    stationlist     = InputDict["stationinfo"]["stationlist"]
    starttime       = InputDict["starttime"]
    endtime         = InputDict["endtime"]
    DL_time_unit    = InputDict["DL_time_unit"]
    DownloadType    = InputDict["DownloadType"]
    fopath          = InputDict["fopath"]

    if mod((endtime - starttime).value,  DL_time_unit) != 0 || (endtime - starttime).value < DL_time_unit
        error("Total download time cannot be devided by Download Time unit; this may cause unexpected result. abort.")
    end

    # calculate start time list (starttimelist) with each Donwload_time_unit
    starttimelist = Utils.get_starttimelist(starttime, endtime, DL_time_unit)
    # generate DLtimestamplist and ststationlist
    DLtimestamplist = Utils.get_timestamplist(starttimelist)

    InputDict["starttimelist"] = starttimelist
    InputDict["DLtimestamplist"] = DLtimestamplist

    #save info into jld2
    jldopen(fopath, "w") do file
        file["info/DLtimestamplist"] = DLtimestamplist;
        file["info/stationlist"] = stationlist;
        file["info/starttime"]   = string(starttime)
        file["info/endtime"]     = string(endtime)
        file["info/DL_time_unit"]= string(DL_time_unit)
    end

    #parallelization by time
    #use all processors
    if DownloadType == "Noise" || DownloadType == "noise"
        S = pmap(x -> seisdownload_NOISE(x, InputDict), 1:length(starttimelist))

    elseif  DownloadType == "Earthquake" || DownloadType == "earthquake"
        println("Download type Earthquake Not implemented.")
    else
        println("Download type is not known (chose Noise or Earthquake).")
    end

    # save data to jld2
    file = jldopen(fopath, "r+")
    unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w")

    for ii = 1:size(S)[1] #loop at each starttime
        for jj = 1:size(S[1])[1] #loop at each station id

            requeststr =S[ii][jj].id
            varname = joinpath(DLtimestamplist[ii], requeststr)
            #save_SeisData2JLD2(fopath, varname, S[ii][jj])
            file[varname] = S[ii][jj]

            if S[ii][jj].misc["dlerror"] == 1
                unavalilablefile[varname] = S[ii][jj]
            end
        end
    end
    JLD2.close(file)
    JLD2.close(unavalilablefile)


    println("Downloading and Saving data is successfully done.\njob ended at "*string(now()))

end

"""
    seisdownload(startid)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `startid`         : start time id in starttimelist
- `InputDict::Dict` : dictionary which contains request information
"""
function seisdownload_NOISE(startid, InputDict::Dict)

    #stationlist
    stationlist     = InputDict["stationinfo"]["stationlist"]
    datacenter      = InputDict["stationinfo"]["stationdatacenter"]
    src             = InputDict["stationinfo"]["stationsrc"]
    starttime       = InputDict["starttime"]
    endtime         = InputDict["endtime"]
    DL_time_unit    = InputDict["DL_time_unit"]
    pre_filt        = InputDict["pre_filt"]

    #make stlist at all processors

    starttimelist = InputDict["starttimelist"]
    timestamplist = InputDict["DLtimestamplist"]

    #show progress
    if starttimelist[startid][end-8:end] == "T00:00:00"
        println("start downloading $(starttimelist[startid])")
    end

    S = SeisData(length(stationlist))
    for i = 1:length(stationlist)
        #---download data---#
        requeststr = stationlist[i]
        argv = [datacenter[i], requeststr, starttimelist[startid], DL_time_unit, 0, src[i], false, "$requeststr.$startid.xml"]
        ex = :(get_data($(argv[1]), $(argv[2]), s=$(argv[3]), t=$(argv[4]), v=$(argv[5]), src=$(argv[6]), w=$(argv[7]), xf=$(argv[8])))
        Stemp = check_and_get_data(ex, requeststr)

        if Stemp.misc[1]["dlerror"] == 0
            Remove_response_obspy.remove_response_obspy!(Stemp, "$requeststr.$startid.xml", pre_filt=pre_filt, output="VEL")
            rm("$requeststr.$startid.xml")
        else
            rm("$requeststr.$startid.xml")
        end

        S[i] = Stemp[1]
    end

    return S
end


"""
    check_and_get_data(ex::Expr, requeststr::String)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `ex::Expr`        : expression of get data includin all request information

# Output
- `S::SeisData`     : downloaded SeisData
- `requeststr::String`     : request channel (e.g. "BP.LCCB..BP1")
"""
function check_and_get_data(ex::Expr, requeststr::String)
   try
       S = eval(ex);
       S.misc[1]["dlerror"] = 0
       return S

   catch
       S = SeisData(1)
       S.misc[1]["dlerror"] = 1
       S.id[1] = requeststr
       note!(S, 1, "station is not available for this request.")
       return S
   end
end

end
