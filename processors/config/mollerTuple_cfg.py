import HpstrConf
import sys
import baseConfig as base
from baseConfig import bfield

base.parser.add_argument("-TS", "--trackstate", type=str, dest="trackstate",
                         help="Specify Track State | 'AtECal', 'AtTarget'. Default is origin ",  metavar="trackstate", default="")

base.parser.add_argument("-r", "--rawHits", type=int, dest="rawHits",
                         help="Keep raw svt hits: 1=yes", metavar="rawHits", default=0)

options = base.parser.parse_args()

# Use the input file to set the output file name
lcio_file = options.inFilename
root_file = options.outFilename

print('LCIO file: %s' % lcio_file)
print('Root file: %s' % root_file)

p = HpstrConf.Process()

# p.max_events = 1000
p.run_mode = 0
p.skip_events = options.skip_events
p.max_events = options.nevents

# Library containing processors
p.add_library("libprocessors")

###############################
#          Processors         #
###############################
header = HpstrConf.Processor('header', 'EventProcessor')
track = HpstrConf.Processor('track', 'TrackingProcessor')
vtx = HpstrConf.Processor('vtx', 'VertexProcessor')
cvtx = HpstrConf.Processor('cvtx', 'VertexProcessor')
mcpart = HpstrConf.Processor('mcpart', 'MCParticleProcessor')
fsp = HpstrConf.Processor('fps', 'FinalStateParticleProcessor')

###############################
#   Processor Configuration   #
###############################
# Event
header.parameters["debug"] = 0
header.parameters["headCollRoot"] = "EventHeader"
header.parameters["trigCollLcio"] = "TriggerBank"
header.parameters["rfCollLcio"] = "RFHits"
header.parameters["vtpCollLcio"] = "VTPBank"
header.parameters["vtpCollRoot"] = "VTPBank"
header.parameters["tsCollLcio"] = "TSBank"
header.parameters["tsCollRoot"] = "TSBank"

# Tracking
track.parameters["debug"] = 0
track.parameters["trkCollLcio"] = 'KalmanFullTracks'
track.parameters["trkCollRoot"] = 'KalmanFullTracks'
track.parameters["kinkRelCollLcio"] = ''
track.parameters["trkRelCollLcio"] = 'KFTrackDataRelations'
track.parameters["trkhitCollRoot"] = 'SiClustersOnTrack'
track.parameters["hitFitsCollLcio"] = 'SVTFittedRawTrackerHits'
track.parameters["rawhitCollRoot"] = 'SVTRawHitsOnTrack_KF'

track.parameters["bfield"] = bfield[str(options.year)]

# Vertex
vtx.parameters["debug"] = 0 
vtx.parameters["vtxCollLcio"] = 'UnconstrainedMollerVertices'
vtx.parameters["vtxCollRoot"] = 'UnconstrainedMollerVertices'
vtx.parameters["partCollRoot"] = 'ParticlesOnMollerVertices'
vtx.parameters["kinkRelCollLcio"] = ''
vtx.parameters["trkRelCollLcio"] = 'KFTrackDataRelations'
vtx.parameters["trkhitCollRoot"] = ''
vtx.parameters["hitFitsCollLcio"] = 'SVTFittedRawTrackerHits'
vtx.parameters["rawhitCollRoot"] = ''
vtx.parameters["trackStateLocation"] = options.trackstate 
if options.trackstate == "":
        vtx.parameters["bfield"] = bfield[str(options.year)]
vtx.parameters["mcPartRelLcio"] = 'SVTTrueHitRelations'

# Vertex
cvtx.parameters["debug"] = 0
cvtx.parameters["vtxCollLcio"] = 'TargetConstrainedMollerVertices'
cvtx.parameters["vtxCollRoot"] = 'TargetConstrainedMollerVertices'
cvtx.parameters["partCollRoot"] = 'ParticlesOnMollerVertices'
cvtx.parameters["kinkRelCollLcio"] = ''
cvtx.parameters["trkRelCollLcio"] = 'KFTrackDataRelations'
cvtx.parameters["trkhitCollRoot"] = ''
cvtx.parameters["hitFitsCollLcio"] = 'SVTFittedRawTrackerHits'
cvtx.parameters["rawhitCollRoot"] = ''
cvtx.parameters["trackStateLocation"] = options.trackstate
if options.trackstate == "":
        cvtx.parameters["bfield"] = bfield[str(options.year)]
cvtx.parameters["mcPartRelLcio"] = 'SVTTrueHitRelations'

# MCParticle
mcpart.parameters["debug"] = 0
mcpart.parameters["mcPartCollLcio"] = 'MCParticle'
mcpart.parameters["mcPartCollRoot"] = 'MCParticle'

#FinalStateParticleProcessor
fsp.parameters["debug"] = 0 
fsp.parameters["fspCollLcio"] = "FinalStateParticles_KF" 
fsp.parameters["fspCollRoot"] = "FinalStateParticles_KF"
fsp.parameters["kinkRelCollLcio"] = ""
fsp.parameters["trkRelCollLcio"] = "KFTrackDataRelations"

if(options.rawHits==1):
    fsp.parameters["trkhitCollRoot"] = "fspOnTrackHits"
    fsp.parameters["rawhitCollRoot"] = "fspOnTrackRawHits"
    fsp.parameters["hitFitsCollLcio"] = "SVTFittedRawTrackerHits"
else:
    fsp.parameters["trkhitCollRoot"] = "fspOnTrackHits"
    fsp.parameters["rawhitCollRoot"] = ""
    fsp.parameters["hitFitsCollLcio"] = ""


sequence = [header, track, vtx, cvtx]

# If MC, get MCParticles
if (not options.isData):
    sequence.append(mcpart)

sequence.append(fsp)

p.sequence = sequence

p.input_files = lcio_file
p.output_files = [root_file]

p.printProcess()

