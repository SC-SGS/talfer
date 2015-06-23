#ifndef GAME_HPP_
#define GAME_HPP_

#include <time.h>
#include <sstream>
#include <string>
#include "SDL.h"
#include "SDL_mixer.h"
#include "CFont.hpp"
#include "CStage_VideoOutput.hpp"
#include "CParticle.hpp"
#include "CParameters.hpp"

class Game {
private:
	int points;
	time_t start_time, end_time;
//	time_t start_pause;
	int elapse_seconds;
//	int paused_seconds;

	int game_length; // in seconds

	/* Mix_Music actually holds the music information.  */

	Mix_Music *music;
	Mix_Chunk *bumped;
	Mix_Chunk *destroyed;
	Mix_Chunk *points_sound;
	Mix_Chunk *spawn_sound;

	CParameters* cParameters;

	CStage_VideoOutput* vo;
	std::string points_s;
	std::string rest_time_s;

public:
	bool global_pause;
	bool has_ended;

public:
	void setCParameters(CParameters* i_cParameters) {
		cParameters = i_cParameters;
	}

	static Game* getInstance() {
		static Game instance;

		return &instance;
	}

public:
	void reset() {
		points = 0;
		points_s = "Points: 0";
		rest_time_s = "";
		global_pause = false;
		has_ended = false;
		time(&start_time);
		end_time = 0;
//		paused_seconds = 0;

		play_music();
		//Mix_RewindMusic();
		//Mix_ResumeMusic();
		// play_music();
	}

	void check_if_paused() {
		time_t now;
		time(&now);
		if ((difftime(now, start_time) >= game_length) && !global_pause) {
			end_game();
			global_pause = true;
		}
		if (cParameters->pause) {
			global_pause = true;
		}
		if (global_pause && !cParameters->pause) {
			reset();
			global_pause = false;
//			paused_seconds += difftime(now,start_pause);
//		}else if(!global_pause && cParameters->pause){
//			time(&start_pause);
		}
//		global_pause = cParameters->pause;
	}

private:
	Game() :
			points(0),
//	paused_seconds(0),
// EDIT: set the length of one game (in seconds), 125 is fine because of the length of the song that is played in the background
			game_length(125), bumped(0), destroyed(0), points_sound(0), spawn_sound(
					0), cParameters(0), points_s("Points: 0"), rest_time_s(""), global_pause(
					false), has_ended(false) {
		// Get game start time
		time(&start_time);
		init_audio();
		load_sounds();
		play_music();
	}

private:
	// Do not implement:
	Game(Game const&);
	void operator=(Game const&);

private:
	void load_sounds() {
		bumped = Mix_LoadWAV("data/music/collision.wav");
		if (!bumped) {
			printf("Mix_LoadWAV(\"collision.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		} else {
			Mix_VolumeChunk(bumped, 32);
		}

		destroyed = Mix_LoadWAV("data/music/destroyed.wav");
		if (!destroyed) {
			printf("Mix_LoadWAV(\"destroyed.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		}

		points_sound = Mix_LoadWAV("data/music/points.wav");
		if (!points_sound) {
			printf("Mix_LoadWAV(\"points.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		}

		spawn_sound = Mix_LoadWAV("data/music/spawn.wav");
		if (!spawn_sound) {
			printf("Mix_LoadWAV(\"spawn.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		} else {
			Mix_VolumeChunk(spawn_sound, 32);
		}

		/* Actually loads up the music */
		music = Mix_LoadMUS("data/music/bg-music.wav");
		if (!music) {
			printf("Mix_LoadMUS(\"bg-music.wav\"): %s\n", Mix_GetError());
			// this might be a critical error...
		}
	}

public:
	void particle_destroyed(float mass, int bumps, float alive_time,
			bool target, int x, int y) {
		// Calc new points to be awarded
		int new_points = mass / (bumps + 1) * alive_time;

		if (target) {
			// Giving points
			points += new_points;

			// TODO Give feedback on amount of points
			if (Mix_PlayChannel(3, destroyed, 0) == -1) {
				printf("Mix_PlayChannel: %s\n", Mix_GetError());
				// well, there's no music, but most games don't break without music...
			}

			// Update point string
			std::stringstream ss;
			ss.flush();
			ss << "Points: " << points;
			points_s = ss.str();
			vo->change_Points(&points_s);

			// Show +Points
			std::stringstream ss2;
			ss2.flush();
			ss2 << "+" << (int) new_points;
			vo->push_Text(ss2.str(), 15, x, y, 20);
		}
	}

	void particle_bumped() {
		if (Mix_PlayChannel(2, bumped, 0) == -1) {
			printf("Mix_PlayChannel: %s\n", Mix_GetError());
			// well, there's no music, but most games don't break without music...
		}
	}

	void speed(float min, float max) {
		// TODO
		// give feedback
		std::cout << "Particle Speeds min:" << min << " max:" << max
				<< std::endl;
	}

	void gesture(std::string /*name*/, float /*direction_x*/,
			float /*direction_y*/) {
		// TODO
		// Displays gesture recognition feedback
	}

	void particle_created(float /*mass*/) {
		if (Mix_PlayChannel(1, spawn_sound, 0) == -1) {
			printf("Mix_PlayChannel: %s\n", Mix_GetError());
			// well, there's no music, but most games don't break without music...
		}
	}

	// particle mass is then calculated as radius^2
	void setNextParticle(CParticle * p) {
		p->color = rand() % 4;
		p->radius = 6 + (rand() % 6);
	}

	/*static void end_game() {
	 Game* g = getInstance();

	 // Get game end time
	 time(&g->end_time);
	 // Calculate elapsed seconds
	 g->elapse_seconds = difftime(g->start_time, g->end_time);

	 g->vo->push_Text("Game Over!", -1, 75, 250, 100);

	 // Halt game
	 g->global_pause = true;

	 g->has_ended=true;
	 }*/

	static void end_game() {
		Game* g = getInstance();

		// Get game end time
		time(&g->end_time);
		// Calculate elapsed seconds
		g->elapse_seconds = difftime(g->start_time, g->end_time);

		// Reset game
//		g->reset();

		// Halt game
		g->cParameters->pause = true;
	}

private:
	void init_audio() {
		/* We're going to be requesting certain things from our audio
		 device, so we set them up beforehand */
		int audio_rate = 22050;
		Uint16 audio_format = AUDIO_S16; /* 16-bit stereo */
		int audio_channels = 2;
		int audio_buffers = 4096;

		/* This is where we open up our audio device.  Mix_OpenAudio takes
		 as its parameters the audio format we'd /like/ to have. */
		if (Mix_OpenAudio(audio_rate, audio_format, audio_channels,
				audio_buffers)) {
			printf("Unable to open audio!\n");
			exit(1);
		}

		// allocate 16 mixing channels
		Mix_AllocateChannels(16);
	}

	void play_music() {
		/* This begins playing the music - the first argument is a
		 pointer to Mix_Music structure, and the second is how many
		 times you want it to loop (use -1 for infinite, and 0 to
		 have it just play once) */
		if (Mix_PlayMusic(music, 0) == -1) {
			printf("Mix_PlayMusic: %s\n", Mix_GetError());
			// well, there's no music, but most games don't break without music...
		}

		//Mix_HookMusicFinished(end_game); // EDIT: commented
	}

public:
	void setVideoOutput(CStage_VideoOutput* vo) {
		this->vo = vo;
		this->updateVideoOutput();
	}

	void update_Time() // updates the rest time string
	{
		if (!global_pause) {
			time_t actual_time;
			if (!end_time)
				time(&actual_time);
			else
				actual_time = end_time;
			int remaining_seconds = game_length
					- difftime(actual_time, start_time);
			int min = remaining_seconds / 60;
			int sec = remaining_seconds - (remaining_seconds / 60) * 60;
			std::stringstream ss;
			ss.flush();
			ss << "Time: ";
			if (min < 10) {
				ss << "0";
			}
			ss << min << ":";
			if (sec < 10) {
				ss << "0";
			}
			ss << sec;
			rest_time_s = ss.str();
		}
	}

	void updateVideoOutput() {
		vo->change_Points(&points_s);
		update_Time();
		vo->change_Time(&rest_time_s);
		//std::cout<< points_s << std::endl;
	}

public:
	void mute_music(bool mute) {
		if (mute) {
			// Turn down volume
			Mix_VolumeMusic(0);
		} else {
			// Turn up volume
			Mix_VolumeMusic(128);
		}
	}
};

#endif
