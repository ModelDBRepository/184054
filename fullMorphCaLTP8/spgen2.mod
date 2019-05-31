TITLE Modified SpikeGenerator of M. Hines
NEURON{ 
  POINT_PROCESS SpGen2
  RANGE y, t01,t02
  RANGE fast_freq, slow_freq, APinburst, startbursting, endbursting
  RANGE noise
}

PARAMETER {
      t01=20 (ms)
      t02=420000 (ms)
	fast_freq	= 100 (/s)	: time between spikes in a burst (msec)
	slow_freq	= 0.00000001 (/s) :0.00000000001:5 (/s)	: burst period (msec)
	APinburst	= 4: 5		: burst length (# spikes)
	startbursting = 30 (ms)	: start of first interburst interval
	endbursting		= 3200 (ms)	: time to stop bursting
	noise		= 0		: amount of randomeaness (0.0 - 1.0)
}

ASSIGNED {
      fast_invl (ms)
      slow_invl (ms)
	y
	burst
	event (ms)
	burst_off (ms)
	burst_on (ms)
	on
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0
	y = -70
    endbursting = endbursting -1e-6 
      if (slow_freq != 0.) {slow_invl= (1/slow_freq )*(1000)}
 if (fast_freq != 0.) {fast_invl= (1/fast_freq )*(1000) }
	burst = 1
       if (t01>0) {
      event= t01
      net_send (event,4)}
 if (t02>0) {
      event= t02
      net_send (event,4)}
      if (startbursting>=0 && APinburst >0 ) {
     event = startbursting
	net_send(event, 3)
	}
}	

FUNCTION interval(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) 	
	}
	if (noise == 0) {
		interval = mean
	}else{
		interval = (1. - noise)*mean + noise*mean*exprand(1)
	}
}

PROCEDURE event_time() {
	if (slow_freq == 0 || (burst != 0. && APinburst > 1)) {
		event = event + interval(fast_invl)
if (slow_freq > 0 && event > burst_on + burst_off) {
			burst = 0.
		}
	}else{
		burst = 1.
: if slow_invl from beginning of burst to beginning of burst
		event = event + interval(slow_invl - (APinburst-1)*fast_invl)
: use following if slow_invl is interburst interval
		:event = event + interval(slow_invl)
		burst_on = event
		burst_off = interval((APinburst - 1)*fast_invl)-1e-6
	}
	if (event > endbursting) {
		on = 0
	}
}

NET_RECEIVE (w) {
:printf("Pregen receive t=%g flag=%g\n", t, flag) 
       if (flag == 3){
      if (on == 0){
       on =1
        event = t
 burst_on = event
		burst_off = interval((APinburst - 1)*fast_invl)-1e-6
         net_send (0,1)
       } }
      if (flag == 4)
{
     if (on == 0) {
      on = 1
      event = t
       net_send (0,1)}}
	if (flag == 1 && on == 1) {
		y = 20
		net_event(t)
       if (t>= startbursting && t<= endbursting)
		{
            event_time()
		net_send(event - t, 1)
           } else {on=0 }
		net_send(.1, 2)
	}
	if (flag == 2) {
		y = -70
	}
}


