#include "zmqextras.h"
#include <zmq.h>
#include <stdlib.h>
#ifdef __WINDOWS__
#include <process.h>
#else
#include <pthread.h>
#endif

struct czmq_thread_ctx {
  czmq_thread_proc *proc;
  void *arg;
  void *ctx;
  void *pipe;
};

#if defined (__WINDOWS__)
static unsigned __stdcall czmq_thread_entry_point(void *arg)
{
  czmq_thread_ctx *thctx = (czmq_thread_ctx *) arg;
  thctx->proc(thctx->arg, thctx->ctx, thctx->pipe);
  free(thctx);
  _endthreadex(0);
  return 0;
}
#else
static void *czmq_thread_entry_point(void *arg)
{
  czmq_thread_ctx *thctx = (czmq_thread_ctx *) arg;
  thctx->proc(thctx->arg, thctx->ctx, thctx->pipe);
  free(thctx);
  return 0;
}
#endif

void *czmq_thread_fork(void *ctx, czmq_thread_proc *proc, void *arg)
{
  char endpoint[256];
  void *pipe = zmq_socket(ctx, ZMQ_PAIR);
  snprintf(endpoint, sizeof(endpoint), "inproc://zctx-pipe-%p", pipe);
  zmq_bind(pipe, endpoint);
  
  czmq_thread_ctx *thctx = (czmq_thread_ctx*)malloc(sizeof(czmq_thread_ctx));
  thctx->ctx = ctx;
  thctx->arg = arg;
  thctx->proc = proc;
  thctx->pipe = zmq_socket(ctx, ZMQ_PAIR);
  zmq_connect(thctx->pipe, endpoint);
  
#if defined (__WINDOWS__)
  HANDLE hThread = (HANDLE)_beginthreadex (NULL, 0, &czmq_thread_entry_point, thctx, CREATE_SUSPENDED, NULL);
  int priority = GetThreadPriority(GetCurrentThread());
  SetThreadPriority(hThread, priority);
  ResumeThread (hThread);
  CloseHandle (hThread);
#else
  pthread_t thread;
  pthread_create(&thread, NULL, czmq_thread_entry_point, thctx);
  pthread_detach(thread);    
#endif
 
  return pipe;
}

bool czmq_poll(void *self, int msecs)
{
  zmq_pollitem_t items [] = { { self, 0, ZMQ_POLLIN, 0 } };
  int rc = zmq_poll (items, 1, msecs);
  return rc != -1 &&
         (items [0].revents & ZMQ_POLLIN) != 0;
}

void czmq_signal(void *self)
{
  zmq_send(self, "", 0, 0);
}

void czmq_wait(void *self)
{
  zmq_msg_t msg;
  zmq_msg_init(&msg);
  zmq_recvmsg(self, &msg, 0);  
  zmq_msg_close(&msg);    
}
