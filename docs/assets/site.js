// This script runs early (not deferred) to set scroll restoration
if ('scrollRestoration' in history) {
  history.scrollRestoration = 'manual';
}

// Ensure we are at the top as early as possible
window.scrollTo(0, 0);

// Reset scroll position to top before leaving the page
// This ensures that if the browser tries to "restore" the scroll position 
// on the next page load (even in manual mode), it restores "0".
window.addEventListener('beforeunload', () => {
  window.scrollTo(0, 0);
});

document.addEventListener('DOMContentLoaded', () => {
  const toggle = document.querySelector(".nav-toggle");
  const sidebar = document.getElementById("sidebar");

  if (toggle && sidebar) {
    toggle.addEventListener("click", () => {
      const open = sidebar.classList.toggle("is-open");
      toggle.setAttribute("aria-expanded", open ? "true" : "false");
    });
  }

  // Handle links to top or internal anchors
  document.addEventListener('click', (e) => {
    const link = e.target.closest('a');
    if (!link) return;
    const href = link.getAttribute('href');
    if (!href) return;

    if (href === '#' || href === '#top' || href === '#site-top') {
      e.preventDefault();
      window.scrollTo({ top: 0, left: 0, behavior: 'smooth' });
      history.pushState(null, null, ' ');
    }
  });

  // Final fallback for page load
  window.addEventListener('load', () => {
    if (!window.location.hash || window.location.hash === '#' || window.location.hash === '#top') {
      window.scrollTo(0, 0);
    }
  });
});
