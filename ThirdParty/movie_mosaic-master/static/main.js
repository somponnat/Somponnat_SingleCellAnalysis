jQuery(document).ready(
  function ($) {

    /* based on http://stackoverflow.com/a/5657492/2674930 */
    var Set = function Set() {
      var list = [];

      var contains;
      this.contains = contains = function(x) {
        return list.indexOf(x) >= 0;
      }

      var put;
      this.put = put = function(x) {
        if (!contains(x)) {
          list.push(x);
        }
        return this;
      }

      var remove;
      this.remove = remove = function(x) {
        var idx = list.indexOf(x);
        if (idx >= 0) {
          list.splice(idx,1);
        }
        return this;
      }

      var first;
      this.first = first = function() {
        return list.length === 0 ? undefined : list[0];
      }

      return this;
    }

    // Track currently-playing videos for sync code below.
    var active_videos = Set();

    var $cells = $('.lookup-table td');

    $cells.click(function() {
      // Find the td's identifying class (starts with "cell-") and use it
      // to look up the corresponding popup element.
      var cell_class = $.grep(this.className.split(/\s+/), function(c) {
        return c.indexOf('cell-') === 0;
      })[0];
      var $popup = $('#popup-' + cell_class);

      $popup.dialog({
        resizable: false,
        show: 'blind',
        width: 'auto',
        height: 'auto',
        position:  {my: 'left top', at: 'left top ', of: window},
        open: function (event) {
          /*
          Implement semi-semantic lazy loading of media elements.
          =====

          This page has many initially hidden jquery-ui dialogs
          containing <img> and <video> tags, but we want the browser
          to put off fetching the media content until it's actually
          displayed. The typical img trick of converting it to a css
          background is undesirable as the image is the actual
          content. Plus that wouldn't help with the videos.

          Here we use an idiom in which we include the media tag as
          usual except for stripping off the src attribute. We then
          put an <a> tag with a class of "media-ondemand"
          immediately preceeding the media tag, and set its href
          attribute to the value of the stripped src attribute from
          the media tag. Here on dialog open we reverse that
          transformation, leaving only fully functional media
          elements behind.
           */
          $('a.media-ondemand', event.target).each(function () {
            var $link = $(this);
            var $media = $link.next();
            // Copy the <a>'s href to the media's src.
            $media.attr('src', $link.attr('href'));
            // Delete the <a> entirely.
            $link.remove();
            // For video tags, create a mediaelement.js player.
            if ($media.prop('nodeName') === 'VIDEO') {
              $media.mediaelementplayer({
                // Ensure the Flash fallback does antialiasing.
                enablePluginSmoothing: true,
                // Allow multiple videos playing at once.
                pauseOtherPlayers: false,
              });
            }
          });
        },
        close: function (event) {
          // On dialog close, pause and rewind all videos.
          $('video', event.target).each(function () {
            try {
              this.pause();
              this.player.setCurrentTime(0);
            } catch (e) {
              // Silence exceptions that happen here if video file
              // failed to load.
            }
          });
        },
      });

    });

    $('video', $cells).bind({
      // Handlers to keep videos synced
      play: function() {
        sync_to_master(this);
        active_videos.put(this);
      },
      pause: function() {
        active_videos.remove(this);
      },
      timeupdate: function() {
        if (this.currentTime === 0 && active_videos.contains(this)) {
          sync_to_master(this);
        }
      },
    });

    function sync_to_master(slave) {
      var master = active_videos.first();
      if (master !== slave && master !== undefined) {
        slave.currentTime = master.currentTime;
      }
    }

  }
);

