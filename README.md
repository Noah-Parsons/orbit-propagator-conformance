# P26 — Orbit propagator conformance

**Field:** Astrodynamics  
**Size:** Large (a main project)  
**Time:** 2-3 months

## Question
Do open-source orbit propagators (poliastro/hapsira, Orekit, Basilisk) agree with analytical Kepler and J2 secular rates, and do they fail silently?

## Why it matters
It carries your conformance method into a new field where errors cost spacecraft.

## Builds on
dynamics-conformance architecture

## Method
1. Derive closed-form references (Kepler, J2 secular drift, Lambert test cases)
2. Adapters for each propagator
3. Pre-register tolerances; run; publish a page like dynamics-conformance

## Deliverable
Second conformance service + report.

## Risks
Frame and time-system conventions cause false alarms, so use your constant-discrepancy diagnostic.

## Before you start: novelty check
Search Google Scholar / arXiv / GitHub for: _orbit propagator verification benchmark open source_. If the question is already answered, narrow it or pick another project.

## Possible outlet
AIAA/AAS student conference; JOSS

## Log
_(dates, decisions, results)_
