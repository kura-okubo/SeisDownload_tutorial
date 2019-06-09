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
- `MAX_MEM_PER_CPU` : maximum available memory for 1 cpu [GB] (default = 1.0GB)
"""
function ParallelSeisrequest(NP::Int, InputDict::Dict; MAX_MEM_PER_CPU::Float64=1.0)

    Utils.initlogo()

	DownloadType    = InputDict["DownloadType"]

    if DownloadType == "Noise" || DownloadType == "noise"

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

        #Test download to evaluate use of memory and estimate download time.
		println("-------TEST DONWLOAD START-----------")
        t = @elapsed Stest = seisdownload_NOISE(1, InputDict) #[s]
        mem_per_requestid = 1.2 * sizeof(Stest) / 1073741824.0 #[GB] *for the safty, required memory is multiplied by 1.2

        max_num_of_processes_per_parallelcycle = floor(Int64, MAX_MEM_PER_CPU/mem_per_requestid)
        estimated_downloadtime = now() + Second(round(10 * t * length(starttimelist) / NP))

		println(mem_per_requestid)
		println(max_num_of_processes_per_parallelcycle)
		println("-------DOWNLOAD STATS SUMMARY------")

		println(@sprintf("Number of processors is %d.", NP))
		println(@sprintf("Total download size will be %4.8f [GB].", mem_per_requestid * length(starttimelist)))
		println(@sprintf("Download will finish at %s.", round(estimated_downloadtime, Dates.Second(1))))

		println("-------START DOWNLOADING-----------")

        if max_num_of_processes_per_parallelcycle < 1
            error("Memory allocation is not enought (currently $MAX_MEM_PER_CPU [GB]). Please inclease MAX_MEM_PER_CPU or decrease number of stations")
        end

        if max_num_of_processes_per_parallelcycle >= length(starttimelist)

            S = pmap(x -> seisdownload_NOISE(x, InputDict), 1:length(starttimelist))

            # save data to jld2
            file = jldopen(fopath, "r+")
            unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")

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

        else

            #parallelization by time
    	    pitr = 1

    	    # progress bar

    	    while pitr <=  length(starttimelist)

    	        startid1 = pitr
    	        startid2 = pitr + max_num_of_processes_per_parallelcycle - 1

    	        if startid2 <= length(starttimelist)
    	            #use all processors
	                S = pmap(x -> seisdownload_NOISE(x, InputDict), startid1:startid2)

    	        else
    	            #use part of processors
	                startid2 = startid1 + mod(length(starttimelist), NP) - 1
	                S = pmap(x -> seisdownload_NOISE(x, InputDict), startid1:startid2)
    	        end

                # save data to jld2
                file = jldopen(fopath, "r+")
                unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")

                for ii = 1:size(S)[1] #loop at each starttime
                    for jj = 1:size(S[1])[1] #loop at each station id

                        requeststr =S[ii][jj].id
                        varname = joinpath(DLtimestamplist[startid1+ii-1], requeststr)
                        #save_SeisData2JLD2(fopath, varname, S[ii][jj])
                        file[varname] = S[ii][jj]

                        if S[ii][jj].misc["dlerror"] == 1
                            unavalilablefile[varname] = S[ii][jj]
                        end
                    end
                end

                JLD2.close(file)
                JLD2.close(unavalilablefile)


    	        pitr += max_num_of_processes_per_parallelcycle

                println("pitr: $pitr")
            end
        end

    elseif  DownloadType == "Earthquake" || DownloadType == "earthquake"

		method		    = InputDict["method"]
		event		    = InputDict["event"]
		reg			    = InputDict["reg"]
		fopath          = InputDict["fopath"]

		#save info into jld2
		jldopen(fopath, "w") do file
			file["info/method"]  = method;
			file["info/event"]   = event;
			file["info/reg"]     = reg
			file["info/fopath"]  = fopath
		end

		#Test download to evaluate use of memory and estimate download time.
		println("-------TEST DONWLOAD START-----------")
		t = @elapsed Stest = seisdownload_EARTHQUAKE(1, InputDict) #[s]
		mem_per_requestid = 1.2 * sizeof(Stest) / 1073741824.0 #[GB] *for the safty, required memory is multiplied by 1.2

		max_num_of_processes_per_parallelcycle = floor(Int64, MAX_MEM_PER_CPU/mem_per_requestid)
		estimated_downloadtime = now() + Second(round(10 * t * length(event) / NP))

		println(mem_per_requestid)
		println(max_num_of_processes_per_parallelcycle)
		println("-------DOWNLOAD STATS SUMMARY------")

		println(@sprintf("Number of processors is %d.", NP))
		println(@sprintf("Total download size will be %4.8f [GB].", mem_per_requestid * length(event)))
		println(@sprintf("Download will finish at %s.", round(estimated_downloadtime, Dates.Second(1))))

		println("-------START DOWNLOADING-----------")

		if max_num_of_processes_per_parallelcycle < 1
			error("Memory allocation is not enought (currently $MAX_MEM_PER_CPU [GB]). Please inclease MAX_MEM_PER_CPU or decrease number of stations")
		end

		if max_num_of_processes_per_parallelcycle >= length(event)

			S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), 1:length(event))

			# save data to jld2
			file = jldopen(fopath, "r+")
			unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")

			for ii = 1:size(S)[1] #loop at each starttime
				varname = joinpath("event",InputDict["event"][ii]["origin"]["time"][1:end-4])
				file[varname] = S[ii]
			end

			JLD2.close(file)
			JLD2.close(unavalilablefile)

		else

			#parallelization by time
			pitr = 1

			# progress bar

			while pitr <=  length(event)

				startid1 = pitr
				startid2 = pitr + max_num_of_processes_per_parallelcycle - 1

				if startid2 <= length(event)
					#use all processors
					S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid2)

				else
					#use part of processors
					startid2 = startid1 + mod(length(event), NP) - 1
					S = pmap(x -> seisdownload_EARTHQUAKE(x, InputDict), startid1:startid2)
				end

				# save data to jld2
				file = jldopen(fopath, "r+")
				unavalilablefile = jldopen(join([fopath[1:end-5], "_unavailablestations.jld2"]), "w+")


				for ii = 1:size(S)[1] #loop at each starttime
					varname = joinpath("event",InputDict["event"][startid1+ii-1]["origin"]["time"][1:end-4])
					file[varname] = S[ii]
				end

				JLD2.close(file)
				JLD2.close(unavalilablefile)

				pitr += max_num_of_processes_per_parallelcycle

				println("pitr: $pitr")
			end
		end


    else
        println("Download type is not known (chose Noise or Earthquake).")
    end

    println("Downloading and Saving data is successfully done.\njob ended at "*string(now()))
    return 0

end

"""
    seisdownload_NOISE(startid, InputDict::Dict)

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
    seisdownload_EARTHQUAKE(startid, InputDict::Dict)

Download seismic data, removing instrumental response and saving into JLD2 file.

# Arguments
- `startid`         : start time id in starttimelist
- `InputDict::Dict` : dictionary which contains request information
"""
function seisdownload_EARTHQUAKE(startid, InputDict::Dict)

	method		    = InputDict["method"]
	event		    = InputDict["event"]
	reg			    = InputDict["reg"]
    pre_filt        = InputDict["pre_filt"]

    #show progress
    if mod(startid, round(0.1*length(event))+1) == 0
        println("start downloading event number: $startid")
    end
    S = SeisData()

    #---download data---#
	for j = 1:length(event[startid]["pickphase"])
		net = event[startid]["pickphase"][j]["net"]
		sta = event[startid]["pickphase"][j]["sta"]
		loc = event[startid]["pickphase"][j]["loc"]
		cha = event[startid]["pickphase"][j]["cha"]

		starttime = event[startid]["pickphase"][j]["starttime"]
		endtime = event[startid]["pickphase"][j]["endtime"]
    	requeststr = join([net,sta,loc,cha], ".")

		if !isempty(InputDict["reg"])
			# request with lat-lon box
    		argv = [method, requeststr, starttime, endtime, InputDict["reg"], 0, false, "$requeststr.$startid.xml"]
		    ex = :(get_data($(argv[1]), $(argv[2]), s=$(argv[3]), t=$(argv[4]), reg=$(argv[5]), v=$(argv[6]), w=$(argv[7]), xf=$(argv[8])))
		    Stemp = check_and_get_data(ex, requeststr)
		else
			# request with lat-lon box
    		argv = [method, requeststr, starttime, endtime, 0, false, "$requeststr.$startid.xml"]
		    ex = :(get_data($(argv[1]), $(argv[2]), s=$(argv[3]), t=$(argv[4]), v=$(argv[5]), w=$(argv[6]), xf=$(argv[7])))
		    Stemp = check_and_get_data(ex, requeststr)
		end

	    if Stemp.misc[1]["dlerror"] == 0
	        Remove_response_obspy.remove_response_obspy!(Stemp, "$requeststr.$startid.xml", pre_filt=pre_filt, output="VEL")
	        rm("$requeststr.$startid.xml")
	    else
	        rm("$requeststr.$startid.xml")
	    end

		append!(S, Stemp)
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
