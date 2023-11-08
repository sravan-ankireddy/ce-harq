# ce-harq
MATLAB code for the ISIT'23 paper: [Compressed Error HARQ: Feedback Communication on Noise-Asymmetric Channels](https://ieeexplore.ieee.org/abstract/document/10206766)

Abstract: In modern communication systems with feedback, there are increasingly more scenarios where the transmitter has much less power than the receiver (e.g., medical implant devices), which we refer to as noise-asymmetric channels. 
For such channels, the feedback link is of higher quality than the forward link. However, feedback schemes for cellular communications, such as hybrid ARQ, do not fully utilize the high-quality feedback link. 
To this end, we introduce Compressed Error Hybrid ARQ, a generalization of hybrid ARQ tailored for noise-asymmetric channels; the receiver sends its estimated message to the transmitter, and the transmitter harmoniously switches between hybrid ARQ and compressed error retransmission. 
We show that our proposed method significantly improves reliability, latency, and spectral efficiency compared to the conventional hybrid ARQ in various practical scenarios where the transmitter is resource-constrained.
