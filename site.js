/* Shared site behavior: mobile navigation, header state, footer year,
   accessible tabs and floating social links. */
(function () {
  'use strict';

  var root = document.documentElement;
  root.classList.add('js');

  /* Mobile navigation */
  var toggle = document.querySelector('.nav-toggle');
  var menu = document.getElementById('primary-menu');

  if (toggle && menu) {
    var setOpen = function (open) {
      toggle.setAttribute('aria-expanded', String(open));
      menu.classList.toggle('is-open', open);
    };

    toggle.addEventListener('click', function () {
      setOpen(toggle.getAttribute('aria-expanded') !== 'true');
    });

    // Closes the menu when a link inside it is chosen
    menu.addEventListener('click', function (event) {
      if (event.target.closest('a')) setOpen(false);
    });

    // Escape closes the menu and returns focus to the toggle
    document.addEventListener('keydown', function (event) {
      if (event.key === 'Escape' && toggle.getAttribute('aria-expanded') === 'true') {
        setOpen(false);
        toggle.focus();
      }
    });

    // Closes the menu when a click lands outside the header
    document.addEventListener('click', function (event) {
      if (toggle.getAttribute('aria-expanded') === 'true' && !event.target.closest('.site-header')) {
        setOpen(false);
      }
    });

    // Resets the menu when the viewport grows past the mobile breakpoint
    var desktop = window.matchMedia('(min-width: 769px)');
    var onChange = function (mq) { if (mq.matches) setOpen(false); };
    if (desktop.addEventListener) desktop.addEventListener('change', onChange);
    else if (desktop.addListener) desktop.addListener(onChange);
  }

  /* Header shadow once the page scrolls */
  var header = document.querySelector('.site-header');
  if (header) {
    var onScroll = function () {
      header.classList.toggle('is-scrolled', window.scrollY > 8);
    };
    onScroll();
    window.addEventListener('scroll', onScroll, { passive: true });
  }

  /* Current year in the footer */
  var year = String(new Date().getFullYear());
  document.querySelectorAll('[data-current-year]').forEach(function (el) {
    el.textContent = year;
  });

  /* Floating social links stay above the footer */
  var social = document.querySelector('.fixed-social-links');
  var footer = document.querySelector('.site-footer');
  if (social && footer) {
    var placeSocial = function () {
      var overlap = window.innerHeight - footer.getBoundingClientRect().top;
      social.style.bottom = overlap > 0 ? (overlap + 24) + 'px' : '';
    };
    placeSocial();
    window.addEventListener('scroll', placeSocial, { passive: true });
    window.addEventListener('resize', placeSocial);
  }

  /* Tabs (WAI-ARIA tabs pattern with automatic activation) */
  document.querySelectorAll('[role="tablist"]').forEach(function (list) {
    var tabs = Array.prototype.slice.call(list.querySelectorAll('[role="tab"]'));

    var select = function (tab, moveFocus) {
      tabs.forEach(function (t) {
        var selected = t === tab;
        t.setAttribute('aria-selected', String(selected));
        t.setAttribute('tabindex', selected ? '0' : '-1');
        t.classList.toggle('active', selected);
        var panel = document.getElementById(t.getAttribute('aria-controls'));
        if (panel) panel.hidden = !selected;
      });
      if (moveFocus) tab.focus();
      if (history.replaceState && tab.dataset.hash) {
        history.replaceState(null, '', '#' + tab.dataset.hash);
      }
    };

    tabs.forEach(function (tab, index) {
      tab.addEventListener('click', function () { select(tab, false); });
      tab.addEventListener('keydown', function (event) {
        var next = null;
        if (event.key === 'ArrowRight') next = tabs[(index + 1) % tabs.length];
        else if (event.key === 'ArrowLeft') next = tabs[(index - 1 + tabs.length) % tabs.length];
        else if (event.key === 'Home') next = tabs[0];
        else if (event.key === 'End') next = tabs[tabs.length - 1];
        if (next) {
          event.preventDefault();
          select(next, true);
        }
      });
    });

    // Opens the tab named in the URL hash, e.g. resume.html#education
    var fromHash = tabs.filter(function (t) {
      return t.dataset.hash && '#' + t.dataset.hash === window.location.hash;
    })[0];
    if (fromHash) select(fromHash, false);
  });
})();
