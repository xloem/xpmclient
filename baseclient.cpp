//============================================================================
// Name        : baseclient.cpp
// Author      : madMAx
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================


#include "baseclient.h"
#include "sysinfo.h"
#include "prime.h"

#include <iostream>

std::string gClientName;
unsigned gClientID;
unsigned gInstanceID;
const unsigned gClientVersion = 11;
const unsigned gClientTarget = 10;

void* gCtx = 0;
static void* gFrontend = 0;
static void* gServer = 0;
static void* gSignals = 0;
static void* gWorkers = 0;

static unsigned gStaleShareChain = 0;
static std::map<unsigned, time_t> timeMap;

std::string gAddr = "";
proto::Block gBlock;
proto::ServerInfo gServerInfo;

static unsigned gNextReqID = 0;
static Timer gRequestTimer;
static unsigned gLatency = 0;
static bool gHeartBeat = true;
static bool gExit = false;

static std::map<unsigned, int> gSharesSent;
static std::map<unsigned, share_t> gShares;
static BaseClient *gClient;


static bool ConnectBitcoin() {
	
	const proto::ServerInfo& sinfo = gServerInfo;
	printf("Connecting to bitcoin: %s:%d ...\n", sinfo.host().c_str(), sinfo.router());

	int linger = 0;
	char endpoint[256];
	snprintf(endpoint, sizeof(endpoint), "tcp://%s:%d", sinfo.host().c_str(), sinfo.router());	
	zmq_close(gServer);
	gServer = zmq_socket(gCtx, ZMQ_DEALER);
	zmq_setsockopt(gSignals, ZMQ_LINGER, &linger, sizeof(int));
	int err = zmq_connect(gServer, endpoint);
	if(err) {
		printf("ERROR: invalid hostname and/or port.\n");
		return false;
	}
	
	return true;
}

static bool ConnectSignals() {
	
	const proto::ServerInfo& sinfo = gServerInfo;
	printf("Connecting to signals: %s:%d ...\n", sinfo.host().c_str(), sinfo.pub());
	
	int linger = 0;
	char endpoint[256];
	snprintf(endpoint, sizeof(endpoint), "tcp://%s:%d", sinfo.host().c_str(), sinfo.pub());
	zmq_close(gSignals);
	gSignals = zmq_socket(gCtx, ZMQ_SUB);
	zmq_setsockopt (gSignals, ZMQ_LINGER, &linger, sizeof(int));
	int err = zmq_connect(gSignals, endpoint);
	if(err){
		printf("ERROR: invalid hostname and/or port.\n");
		return false;
	}
	
	const char one[2] = {1, 0};
	zmq_setsockopt (gSignals, ZMQ_SUBSCRIBE, one, 1);
	
	return true;
	
}


inline void GetNewReqNonce(proto::Request& req) {
	
	uint32_t limbs[8];
	for(int i = 0; i < 8; ++i)
		limbs[i] = rand();
	
	uint32_t tmp = limbs[0];
	for(int i = 1; i < 7; ++i)
		tmp *= limbs[i];
	
	limbs[7] = -tmp;
	req.set_reqnonce(&limbs[0], sizeof(limbs));
	
}


static void RequestWork() {
	
	if(!gServer)
		return;
	
	//printf("Requesting new work...\n");
	
	proto::Request req;
	req.set_type(proto::Request::GETWORK);
	req.set_reqid(++gNextReqID);
	req.set_version(gClientVersion);
	req.set_height(gBlock.height());
	
	GetNewReqNonce(req);
	Send(req, gServer);
	
	gRequestTimer.reset();
	
}


static void HandleNewBlock(const proto::Block& block) {
	
	if(block.height() > gBlock.height()){
		
		//printf("New Block: height=%d\n", block.height());
		
		gBlock = block;
		gClient->NotifyBlock(block);
		
		RequestWork();
		
	}
	
}


static void HandleNewWork(const proto::Work& work) {
	
	float diff = gRequestTimer.diff();
	gLatency = diff * 1000.;
	
	printf("Work received: height=%d diff=%.8g latency=%dms\n",
			work.height(), GetPrimeDifficulty(work.bits()), gLatency);
	
	gClient->TakeWork(work);
	
}


static void SubmitShare(const proto::Share& share) {
	
	proto::Request req;
	req.set_type(proto::Request::SHARE);
	req.set_reqid(++gNextReqID);
	req.set_version(gClientVersion);
	req.set_height(share.height());
	req.mutable_share()->CopyFrom(share);
	
	GetNewReqNonce(req);
	Send(req, gServer);
	
  timeMap[req.reqid()] = time(0);
	gSharesSent[req.reqid()] = share.length();
	
}


static void PrintStats() {
	
	printf("(ST/INV/DUP):");
	for(std::map<unsigned, share_t>::const_iterator iter = gShares.begin(); iter != gShares.end(); ++iter){
		
		const share_t& share = iter->second;
		printf(" %dx %dch(%d/%d/%d)", share.accepted, iter->first, share.stale, share.invalid, share.duplicate);
		
	}
	
	if(!gShares.size())
		printf(" none");
		
	printf("\n");
	
}


static int HandleReply(void *socket) {
	
	proto::Reply rep;
	Receive(rep, socket);
	
	if(rep.has_errstr())
		printf("Message from server: %s\n", rep.errstr().c_str());
	
	if(rep.has_block()){
		
		HandleNewBlock(rep.block());
		
	}
	
	if(rep.type() == proto::Request::GETWORK){
		
		if(rep.has_work()){
			
			HandleNewWork(rep.work());
			
		}
		
	}else if(rep.type() == proto::Request::SHARE){
		timeMap.erase(rep.reqid());
    
		unsigned length = gSharesSent[rep.reqid()];
		gSharesSent.erase(rep.reqid());

		switch(rep.error()){
		case proto::Reply::NONE:
			printf("Share accepted.\n");
			gStaleShareChain = 0;
			gShares[length].accepted++;  
			break;
		case proto::Reply::INVALID:
			printf("Invalid share.\n");
			gShares[length].invalid++;
			break;
		case proto::Reply::STALE:
			printf("Stale share.\n");
			gShares[length].stale++;
			gStaleShareChain++;
			if (gStaleShareChain >= 4) {
				gStaleShareChain = 0;
				gExit = false;
				return -1;
			}
			break;
		case proto::Reply_ErrType_DUPLICATE:
			printf("Duplicate share.\n");
			gShares[length].duplicate++;
			break;
		default: break;
		}
		
	}else if(rep.type() == proto::Request::STATS){
		
		if(rep.error() == proto::Reply::NONE)
			{}//printf("Statistics accepted.\n");
		else
			printf("Statistics rejected.\n");
		
	}else if(rep.type() == proto::Request::PING){
		
		//printf("Received heartbeat.\n");
		
	}
	
	gHeartBeat = true;
	
	return 0;
	
}


static int HandleSignal(void *socket) {
	
	proto::Signal sig;
	ReceivePub(sig, socket);
	
	if(sig.type() == proto::Signal::NEWBLOCK){
		
		HandleNewBlock(sig.block());
		
	}else if(sig.type() == proto::Signal::SHUTDOWN){
		
		printf("Server is shutting down, reconnecting...\n");
		gExit = false;
		return -1;
		
	}
	
	return 0;
	
}


static int HandleWorkers(void *socket) {
	
	proto::Share share;
	bool ret = Receive(share, socket);
	if(!ret)
		return 0;

	SubmitShare(share);
	
	return 0;
	
}


static int HandleTimer() {
	
	if(!gHeartBeat){
		
		printf("Connection lost, reconnecting...\n");
		gExit = false;
		return -1;
		
	}
	
	{
		proto::Request req;
		req.set_type(proto::Request::STATS);
		req.set_reqid(++gNextReqID);
		req.set_version(gClientVersion);
		req.set_height(gBlock.height());
		
		proto::ClientStats& stats = *req.mutable_stats();
		int ngpus = gClient->GetStats(stats);
		stats.set_addr(gAddr);
		stats.set_name(gClientName);
		stats.set_clientid(gClientID);
		stats.set_instanceid(gInstanceID);
		stats.set_version(gClientVersion);
		stats.set_latency(gLatency);
		stats.set_ngpus(ngpus);
		stats.set_height(gBlock.height());
		
		GetNewReqNonce(req);
		Send(req, gServer);
	}
	
	PrintStats();
	
	gHeartBeat = false;
	
	//ReConnectSignals();
	
	return 0;
	
}


static int TimeoutCheckProc() {  
  time_t currentTime = time(0);
  for (auto &t: timeMap) {
    if (currentTime - t.second >= 5) {
      timeMap.clear();
      printf("Connection lost, reconnecting...\n");
      gExit = false;
      return -1;
    }
  }
  
  return 0;
  
}


int main(int argc, char **argv) {
gBlock.set_height(0);
	gClientName = sysinfo::GetClientName();
	gClientID = sysinfo::GetClientID();
	gInstanceID = gClientID * (unsigned)time(0);
	srand(gInstanceID);
	
	std::string frontHost;
	unsigned frontPort;
	
	Configuration* cfg = Configuration::create();
	try{
		cfg->parse("config.txt");
		frontHost = cfg->lookupString("", "server", "localhost");
		frontPort = cfg->lookupInt("", "port", 6666);
		gAddr = cfg->lookupString("", "address", "");
		gClientName = cfg->lookupString("", "name", gClientName.c_str());
	}catch(const ConfigurationException& ex){
		printf("ERROR: %s\n", ex.c_str());
		printf("hit return to exit...\n");
		std::string line;
		std::getline(std::cin, line);
		exit(EXIT_FAILURE);
	}
	
	if(!gClientName.size())
		gClientName = sysinfo::GetClientName();
	
	printf("madPrimeMiner-v%d.%d\n", gClientVersion/10, gClientVersion%10);
	printf("ClientName = '%s'  ClientID = %u  InstanceID = %u\n", gClientName.c_str(), gClientID, gInstanceID);
	printf("Address = '%s'\n", gAddr.c_str());
	
	if(!gAddr.size()){
		printf("ERROR: address not specified in config.txt\n");
		printf("hit return to exit...\n");
		std::string line;
		std::getline(std::cin, line);
		exit(EXIT_FAILURE);
	}

	gCtx = zmq_ctx_new();
	gWorkers = zmq_socket(gCtx, ZMQ_PULL);
	zmq_bind(gWorkers, "inproc://shares");

  gClient = createClient(gCtx);
  
  bool benchmarkOnly = false;
  if (argc >= 2 && (strcmp(argv[1], "-b") == 0 || strcmp(argv[1], "--benchmark") == 0))
    benchmarkOnly = true;
	gExit = !gClient->Initialize(cfg, benchmarkOnly);
	
	while(!gExit){
		printf("Connecting to frontend: %s:%d ...\n", frontHost.c_str(), frontPort);
		
		gBlock.Clear();
		proto::Reply rep;
		gExit = true;
		
		while(gExit){
			int linger = 0;
			char endpoint[256];
			snprintf(endpoint, sizeof(endpoint), "tcp://%s:%d", frontHost.c_str(), frontPort);
			zmq_close(gFrontend);
			gFrontend = zmq_socket(gCtx, ZMQ_DEALER);
			zmq_setsockopt (gFrontend, ZMQ_LINGER, &linger, sizeof(int));                        
			int err = zmq_connect(gFrontend, endpoint);
			if(err){
				printf("ERROR: invalid hostname and/or port.\n");
				exit(EXIT_FAILURE);
			}
			
			proto::Request req;
			req.set_type(proto::Request::CONNECT);
			req.set_reqid(++gNextReqID);
			req.set_version(gClientVersion);
			req.set_height(0);
			
			GetNewReqNonce(req);
			Send(req, gFrontend);
			
			bool ready = czmq_poll(gFrontend, 3*1000);                        
			if(!ready)
				continue;
			
			Receive(rep, gFrontend);
			
			if(rep.error() != proto::Reply::NONE){
				printf("ERROR: %s\n", proto::Reply::ErrType_Name(rep.error()).c_str());
				if(rep.has_errstr())
					printf("Message from server: %s\n", rep.errstr().c_str());
			}
			
			if(!rep.has_sinfo())
				break;
			
			gServerInfo = rep.sinfo();
			
			bool ret = false;
			ret |= !ConnectBitcoin();
			ret |= !ConnectSignals();
			if(ret)
				break;
			
			gExit = false;
			
		}

		char endpoint[256];
		snprintf(endpoint, sizeof(endpoint), "tcp://%s:%d", frontHost.c_str(), frontPort);
		zmq_disconnect(gFrontend, endpoint);
		
		if(gExit)
			break;

       bool loopActive = true;
       time_t timer1sec = time(0);
       time_t timer1min = time(0);                
       zmq_pollitem_t items[] = {
         {gServer, 0, ZMQ_POLLIN, 0},
         {gSignals, 0, ZMQ_POLLIN, 0},
         {gWorkers, 0, ZMQ_POLLIN, 0}
       };
		
		gHeartBeat = true;
		gExit = true;
		
		if(rep.has_block())
			HandleNewBlock(rep.block());
		else
			RequestWork();
		
		gClient->Toggle();

		while (loopActive) {
			int result = zmq_poll(items, sizeof(items)/sizeof(zmq_pollitem_t), 1000);
			if (result == -1)
				break;
         
			if (result > 0) {
				if (items[0].revents & ZMQ_POLLIN)
					loopActive &= (HandleReply(gServer) == 0);
				if (items[1].revents & ZMQ_POLLIN)
					loopActive &= (HandleSignal(gSignals) == 0);
				if (items[2].revents & ZMQ_POLLIN)
					loopActive &= (HandleWorkers(gWorkers) == 0);
			}
			
			// check timers
			time_t currentTime = time(0);
			if (currentTime - timer1sec >= 1) {
				timer1sec = currentTime;
				loopActive &= (TimeoutCheckProc() == 0);
			}
			
			if (currentTime - timer1min >= 60) {
				timer1min = currentTime;
				loopActive &= (HandleTimer() == 0);
			}
		}                
		
		gClient->Toggle();
		zmq_close(gServer);
		zmq_close(gSignals);
		gServer = 0;
		gSignals = 0;
	}
	
	delete gClient;
	zmq_close(gWorkers);
	zmq_close(gFrontend);
	zmq_ctx_shutdown(gCtx);
	
	return EXIT_SUCCESS;
	
}

BaseClient::BaseClient(void *ctx)
{
  mCtx = ctx;
  mBlockPub = zmq_socket(ctx, ZMQ_PUB);
  mWorkPub = zmq_socket(ctx, ZMQ_PUB);
  mStatsPull = zmq_socket(ctx, ZMQ_PULL);
 
  zmq_bind(mBlockPub, "inproc://blocks");
  zmq_bind(mWorkPub, "inproc://work");
  zmq_bind(mStatsPull, "inproc://stats");

  mPaused = true;
  mNumDevices = 0;
  mStatCounter = 0;
}
